"""This module implements a solver for the Einsten-Vlasov
equations in axial symmetry."""

import os, sys, numpy

from dolfin import *
from ufl.algorithms import extract_coefficients

from models import *

class AdaptiveVlasovPoissonSolver:
    "Adaptive solver for the axisymmetric Vlasov-Poisson equations"

    def __init__(self):

        # Create parameter set
        self.parameters = Parameters(discretization = Parameters("discretization"),
                                     output = Parameters("output"))

        # Set up directories
        library_dir  = os.path.dirname(os.path.abspath(__file__))
        geometry_dir = os.path.join(library_dir, "geometries")
        solution_dir = os.path.join("solutions", "vp")

        # Discretization parameters
        self.parameters.discretization.add("mass", 1.0)
        self.parameters.discretization.add("radius", 25)
        self.parameters.discretization.add("maxiter", 1000)
        self.parameters.discretization.add("tolerance", 1e-3)
        self.parameters.discretization.add("num_steps", 10)
        self.parameters.discretization.add("resolution", 32)
        self.parameters.discretization.add("resolution_3d", 8)
        self.parameters.discretization.add("mesh_prefix", "halfcircle")
        self.parameters.discretization.add("adaptive", False)

        # Output parameters
        self.parameters.output.add("library_directory",  library_dir)
        self.parameters.output.add("geometry_directory", geometry_dir)
        self.parameters.output.add("solution_directory", solution_dir)
        self.parameters.output.add("suffix",             "")
        self.parameters.output.add("plot_solution",      True)
        self.parameters.output.add("save_solution",      True)
        self.parameters.output.add("save_solution_3d",   True)
        self.parameters.output.add("save_point_cloud",   True)

        # FIXME: Move this to common function for all solvers
        # Override some parameters when --hires option is given
        if len(sys.argv) > 1 and sys.argv[1] == "--hires":
            print("Running simulation with high resolution")
            self.parameters.discretization.num_steps = 50
            self.parameters.discretization.resolution = 512
            self.parameters.discretization.resolution_3d = 64
            self.parameters.discretization.tolerance = 1e-4
        if len(sys.argv) > 2:
            solution_dir = sys.argv[2]
            self.parameters.output.solution_directory = solution_dir

    def _create_initial_values(self, V, e0, m, r, z):
        "Create initial values for iteration"

        # Parameters
        E0 = Constant(e0)
        M  = Constant(m)
        r0 = -M / E0

        # Create function
        U = Function(V)

        # Project to finite element space
        project(2*E0 / (1 + sqrt(r**2 + z**2) / r0), V=V, function=U)

        return U

    def _create_forms(self, V, U, rho, r):
        "Create forms"

        # Create trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)

        # Create forms
        a = self._create_lhs(u, v, r)
        L = self._create_rhs(rho, v, r)
        F = self._create_lhs(U, v, r) - L

        return a, L, F

    def _create_lhs(self, u, v, r):
        "Create left-hand side variational form"
        a = -inner(grad(u), grad(v))*r*dx
        return a

    def _create_rhs(self, rho, v, r):
        "Create right-hand side variational form"
        L = 4*pi*rho*v*r*dx
        return L

    def _create_residual(self, u, v, rho, r):
        "Create residual form"
        lhs = self._create_lhs(u, v, r)
        rhs = self._create_rhs(rho, v, r)
        R = lhs - rhs
        return R

    def _create_error_indicators(self, U, rho):
        "Create form for error indicators"

        # Extract mesh
        mesh = U.function_space().mesh()

        # Create bubble function space and scaled bubble test function
        B = FunctionSpace(mesh, "Bubble", 3)
        h = CellSize(mesh)
        x = SpatialCoordinate(mesh)
        r = x[0]
        v = h*h*TestFunction(B)
        #v = TestFunction(B)

        # Create form for error indicator
        eT = self._create_residual(U, v, rho, r)

        eT = U.dx(0)*v*dx


        return eT

    def _create_boundary_conditions(self, V, m, R):
        "Create boundary conditions"

        # Note: can't use boundary markers here since those
        # are no longer valid after refining the mesh...

        R = 0.99*R

        bc  = DirichletBC(V, -m / R, "x[0]*x[0] + x[1]*x[1] > %g" % R)
        bc0 = DirichletBC(V, 0.0,    "x[0]*x[0] + x[1]*x[1] > %g" % R)

        return bc, bc0

    def _refine_mesh(self, mesh, U, rho):
        "Refine mesh based on error indicator form"

        print("Refining mesh...")

        # Create error indicator form
        eT = self._create_error_indicators(U, rho)

        #W = FunctionSpace(mesh, "DG", 0)
        #E = Function(W)
        #assemble(eT, tensor=E.vector())
        #plot(E, interactive=True)

        # Assemble error indicators
        ET = assemble(eT)
        indicators = [(numpy.abs(e), i) for i, e in enumerate(ET.array())]
        indicators = list(reversed(sorted(indicators)))

        # Mark cells based on indicators
        markers = CellFunction("bool", mesh)
        markers.set_all(False)
        dorfler_fraction = 0.3
        dorfler_sum = 0.0
        dorfler_max = numpy.sum(e for (e, i) in indicators)
        for e, i in indicators:
            markers.set_value(i, True)
            dorfler_sum += e
            if dorfler_sum > dorfler_fraction*dorfler_max:
                break

        # Refine mesh
        refined_mesh = refine(mesh, markers)

        return refined_mesh

    def _save_data(self, solution_dir, data):
        "Save data (append to file if any)"

        # Create directory if it does not yet exist
        if not os.path.exists(solution_dir):
            os.makedirs(solution_dir)

        # Name of file
        filename = os.path.join(solution_dir, "data.csv")
        info("Appending data to file %s." % filename)

        # Get keys
        keys = sorted(data.keys())

        # Write header if not written before
        if not os.path.isfile(filename):
            f = open(filename, "w")
            f.write(",".join(str(k) for k in keys) + "\n")
            f.close()

        # Append data
        f = open(filename, "a")
        f.write(",".join("%.16g" % data[k] for k in keys) + "\n")
        f.close()

    def _iterate(self, model, mesh, U,
                 e0, m, R, ansatzes, num_steps):
        "Perform one iteration"

        # FIXME: Get rid of very long parameter list above

        # Create spatial coordinates
        x = SpatialCoordinate(mesh)
        r = x[0]
        z = x[1]

        # Create function space and function
        if U is None:
            V = FunctionSpace(mesh, "Lagrange", 1)
            U = self._create_initial_values(V, e0, m, r, z)
        else:
            V = U.function_space()

        # Set density functional
        _rho = model
        C = Constant(1)
        rho = C*_rho

        # Define unscaled and scaled mass
        _mass = 2*2*pi*_rho*r*dx
        mass  = C*_mass

        # Create boundary conditions
        bc, bc0 = self._create_boundary_conditions(V, m, R)

        # Create forms
        a, L, F = self._create_forms(V, U, rho, r)

        # Initialize all ansatzes
        for ansatz in ansatzes:
            ansatz.set_integration_parameters(num_steps)
            ansatz.read_parameters()

        # Assemble left-hand side and apply boundary conditions
        A = assemble(a)
        b = Vector()
        bc.apply(A)

        # Create linear solver
        preconditioners = [pc for pc in krylov_solver_preconditioners()]
        if "amg" in preconditioners:
            solver = LinearSolver("gmres", "amg")
        else:
            warning("Missing AMG preconditioner, using ILU.")
            solver = LinearSolver("gmres", "ilu")

        # Set linear solver parameters
        solver.parameters["relative_tolerance"] = 1e-9

        # Print solver discretization parameters
        self.print_discretization()

        # Reset computation of radius of support
        for ansatz in ansatzes:
            ansatz.set_potential(U)
            ansatz.reset()

        # Rescale right-hand side to preserve mass
        _m = assemble(_mass)
        C.assign(m / _m)
        info("C = %.15g" % float(C))

        # Assemble right-hand side and apply boundary conditions
        assemble(L, tensor=b)
        bc.apply(b)

        # Solve linear system
        solver.solve(A, U.vector(), b)

        # Compute residual
        f = assemble(F)
        bc0.apply(f)
        residual = norm(f)
        info("||F|| = %.3g" % residual)

        return U, rho, residual

    def solve(self, model, solution=None):
        "Compute solution"

        # Extract all ansatzes from model
        ansatzes = [c for c in extract_coefficients(model) if not isinstance(c, Constant)]

        # Get common model parameters (may be different but look only at first)
        e0 = ansatzes[0].parameters.E0

        # Get discretization parameters
        R             = self.parameters.discretization.radius
        m             = self.parameters.discretization.mass
        maxiter       = self.parameters.discretization.maxiter
        tol           = self.parameters.discretization.tolerance
        num_steps     = self.parameters.discretization.num_steps
        N             = self.parameters.discretization.resolution
        resolution_3d = self.parameters.discretization.resolution_3d
        mesh_prefix   = self.parameters.discretization.mesh_prefix
        adaptive      = self.parameters.discretization.adaptive

        # Get output parameters
        geometry_dir     = self.parameters.output.geometry_directory
        solution_dir     = self.parameters.output.solution_directory
        suffix           = self.parameters.output.suffix
        plot_solution    = self.parameters.output.plot_solution
        save_solution    = self.parameters.output.save_solution
        save_solution_3d = self.parameters.output.save_solution_3d
        save_point_cloud = self.parameters.output.save_point_cloud

        # FIXME: May not be needed
        parameters.allow_extrapolation = True

        # Set file names for mesh and markers
        mesh_file = os.path.join(geometry_dir, "%s_%d_%d_mesh.xml.gz" \
                                 % (mesh_prefix, R, N))

        # Initialize mesh and solution
        mesh = Mesh(mesh_file)
        U = None
        rho = None

        # FIXME: Temporary post-processing
        # File for storing adaptive meshes
        mfile = File("solutions/avp/mesh.pvd")
        ufile = File("solutions/avp/U.pvd")
        rfile = File("solutions/avp/rho.pvd")        

        # Main loop
        tic()
        for iter in range(maxiter):

            begin("Iteration %d" % iter)

            # Perform one iteration
            U, rho, residual = self._iterate(model, mesh, U,
                                            e0, m, R, ansatzes, num_steps)

            mfile << U.function_space().mesh()
            ufile << U

            # Refine mesh
            mesh = self._refine_mesh(mesh, U, rho)
            V = FunctionSpace(mesh, "Lagrange", 1)
            U = interpolate(U, V)

            RHO = Function(V)
            project(rho, mesh=mesh, function=RHO)
            rfile << RHO            

            # Check for convergence
            if residual < tol and iter > 0:
                break

        # Print elapsed time
        print(("Iterations finished in %g seconds." % toc()))

        # Check whether iteration converged
        if iter == maxiter - 1:
            error("Iteration did not converge.")
        print()
        print(("Iterations converged to within a tolerance of %g." % tol))
        print(("Number of iterations was %g." % iter))

        # FIXME: Temporary post-processing
        plot(mesh, interactive=True)
        #plot(U, interactive=True)

        # FIXME: Stuff below needs fixing
        return U, None, None

        # Compute some interesting values
        _m = assemble(_mass)
        C.assign(m / _m)
        m  = assemble(mass)
        r0 = max([ansatz.radius_of_support() for ansatz in ansatzes])
        R0   = r0*(1.0 + m/(2.0*r0))**2

        # Store values
        self.data = {"mass": m,
                     "unscaled_mass": _m,
                     "radius_of_support": r0,
                     "areal_radius_of_support": R0}

        # Print data
        self.print_ansatzes(ansatzes)
        self.print_data()

        # Save solutions to file
        if save_solution:

            # Save data
            self._save_data(solution_dir, self.data)

            # Save to XML format for later reuse
            File(os.path.join(solution_dir, "mesh_%d%s.xml" % (R, suffix))) << mesh
            File(os.path.join(solution_dir, "U_%d%s.xml" % (R, suffix)))    << U
            File(os.path.join(solution_dir, "RHO_%d%s.xml" % (R, suffix)))  << RHO

            # Save to PVD format for plotting
            File(os.path.join(solution_dir, "U_%d%s.pvd" % (R, suffix))) << U
            File(os.path.join(solution_dir, "RHO_%d%s.pvd" % (R, suffix))) << RHO

            # Save Schwarzschild solution on an exterior annulus
            #mesh_file = os.path.join(geometry_dir, "halfannulus_%d_mesh.xml.gz" % R)
            #annulus_mesh = Mesh(mesh_file)
            #Z = FunctionSpace(annulus_mesh, "Lagrange", 1)
            #_U_R = interpolate(U_R, Z)
            #File(os.path.join(solution_dir, "U_R_%d%s.pvd" % (R, suffix))) << _U_R

        # Interpolate to create 3D representation of density
        if save_solution_3d:
            print("Computing 3D representation of density")
            n = resolution_3d
            box = BoxMesh(Point(-R, -R, -R), Point(R, R, R), n, n, n)
            V3D = FunctionSpace(box, "Lagrange", 1)
            rho3d = interpolate(Density3D(RHO), V3D)
            File(os.path.join(solution_dir, "RHO3D_%d%s.pvd" % (R, suffix))) << rho3d

        # Create point cloud representation of density
        if save_point_cloud:
            print("Computing point cloud representation of density")
            rho = PointCloud(RHO, R, m, 64, 50000)
            filename = os.path.join(solution_dir, "point_cloud_%d%s.xdmf" % (R, suffix))
            rho.save_data(filename)

        # Plot solutions
        if plot_solution:

            # Plot metric (potential)
            plot(U, title="U")

            # Plot density distribution
            project(density, mesh=mesh, function=RHO)
            plot(RHO, title="Density")

            interactive()

        info("Solver finished")

        return U, RHO, self.data

    def print_ansatzes(self, ansatzes):

        for ansatz in ansatzes:
            print()
            info(ansatz.parameters, True)

    def print_discretization(self):

        print()
        print("Discretization parameters:")
        info(self.parameters.discretization, True)
        print()


    def print_data(self):
        "Pretty-print data from computation"

        r0  = self.data["radius_of_support"]
        R0 = self.data["areal_radius_of_support"]
        M  = self.data["mass"]
        _M = self.data["unscaled_mass"]
        Q  = 2*M / R0

        print()
        print(("  r0  = %.16g \t(Radius of support)"               % r0))
        print(("  R0  = %.16g \t(Radius of support areal coords)"  % R0))
        print(("  Q  = %.16g \t(2*M / R0)"                         % Q))
        print(("  M  = %.16g \t(Total mass)"                       % M))
        print(("  M' = %.16g \t(Unscaled mass)"                    % _M))
        print()

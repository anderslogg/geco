"""This module implements a solver for the Einsten-Vlasov equations in
axial symmetry with or without net angular momentum."""

from solverbase import *

def _lhs(nu, bb, mu, ww, v0, v1, v2, v3, r,
         P00, P11, P33, P03):

    a0 = dot(grad(nu), grad(v0))*r*dx
    a1 = dot(grad(bb), grad(v1))*r*dx - bb.dx(0)*v1*dx + 8*pi*bb*P11*v1*r*dx
    a2 = dot(grad(mu), grad(v2))*r*dx + mu.dx(0)*v2*dx
    a3 = dot(grad(ww), grad(v3))*r*dx - 2.0*ww.dx(0)*v3*dx

    return a0, a1, a2, a3

def _rhs(nu, bb, mu, ww, v0, v1, v2, v3, r,
         P00, P11, P33, P03):

    L0 = - 4*pi*(P00 + P11)*v0*r*dx \
         - 4*pi*(1.0 + (r*bb)**2*exp(-4*nu)*ww**2)*P33*v0*r*dx \
         - 8*pi*exp(-4*nu)*ww*P03*v0*r*dx \
         + (1/bb)*dot(grad(bb), grad(nu))*v0*r*dx \
         - 0.5*exp(-4.0*nu)*(r*bb)**2*dot(grad(ww), grad(ww))*v0*r*dx

    L1 = Constant(0.0)*v1*dx

    L2 = 4*pi*(P00 + P11)*v2*r*dx \
         + 4*pi*((r*bb)**2*exp(-4.0*nu)*ww**2 - 1.0)*P33*v2*r*dx \
         + 8*pi*exp(-4.0*nu)*ww*P03*v2*r*dx \
         - (1/bb)*dot(grad(bb), grad(nu))*v2*r*dx \
         - nu.dx(0)*v2*dx \
         + dot(grad(nu), grad(nu))*v2*r*dx \
         - 0.25*exp(-4.0*nu)*(r*bb)**2*dot(grad(ww), grad(ww))*v2*r*dx

    L3 = - 16*pi/(r*bb)**2*(P03 + ww*(r*bb)**2*P33)*v3*r*dx \
         + 3.0*(1/bb)*dot(grad(bb), grad(ww))*v3*r*dx \
         - 4.0*dot(grad(nu), grad(ww))*v3*r*dx

    return L0, L1, L2, L3

def _flat(m, J):

    parameters = {"m": m, "J": J}

    nu = "-m / sqrt(x[0]*x[0] + x[1]*x[1])"
    bb = "1.0 - m*m / (4.0*(x[0]*x[0] + x[1]*x[1]))"
    mu = "m / sqrt(x[0]*x[0] + x[1]*x[1])"
    ww = "2.0*J / pow(x[0]*x[0] + x[1]*x[1], 1.5)"

    return [Expression(v, degree=1, **parameters) for v in (nu, bb, mu, ww)]

def _init(e0, V):

    parameters = {"E": e0, "r0": 4.0, "c": 0.25, "d": 0.25}

    nu = "-0.5*E / (1.0 + sqrt(x[0]*x[0] + x[1]*x[1]) / r0)"
    bb = "1.0 - c*exp(-d*sqrt(x[0]*x[0] + x[1]*x[1]))"
    mu = "E / (1.0 + sqrt(x[0]*x[0] + x[1]*x[1]) / r0)"
    ww = "0.0"

    return [project(Expression(v, degree=1, **parameters), V) for v in (nu, bb, mu, ww)]

class EinsteinVlasovSolver(SolverBase):
    "Solver for the axisymmetric Einstein-Vlasov equations"

    def __init__(self):
        SolverBase.__init__(self, "ev")

        # Add special parameter
        self.parameters.discretization.add("ang_mom", 0.0)

    def solve(self, model, solution=None):
        "Compute solution"

        # Initialize solve
        self._init_solve()
        
        # Extract all ansatzes from model
        ansatzes = [c for c in extract_coefficients(model) if not isinstance(c, Constant)]

        # Get common model parameters (use first)
        e0 = ansatzes[0].parameters.E0

        # Get discretization parameters
        m         = self.parameters.discretization.mass
        J         = self.parameters.discretization.ang_mom
        maxiter   = self.parameters.discretization.maxiter
        theta     = self.parameters.discretization.theta
        tol       = self.parameters.discretization.tolerance
        num_steps = self.parameters.discretization.num_steps
        R         = self.parameters.discretization.radius
        degree    = self.parameters.discretization.degree
        depth     = self.parameters.discretization.anderson_depth

        # Get output parameters
        plot_iteration = self.parameters.output.plot_iteration

        # Workaround for geometric round-off errors
        parameters.allow_extrapolation = True

        # Generate mesh and create function space
        mesh = self._generate_mesh()
        V = FunctionSpace(mesh, "Lagrange", degree)

        # Create initial data
        U = _init(e0, V)
        RHO = Function(V)

        # Override initial data if provided
        if solution is not None:
            info("Reusing function space and initial data.")
            U = solution[0:4]
            V = U[0].function_space()
            RHO = Function(V)
            mesh = V.mesh()

        # Create spatial coordinates
        x = SpatialCoordinate(mesh)
        r = x[0]
        z = x[1]

        # Extract solution components
        NU, BB, MU, WW = U

        # Create asymptotically flat solutions for boundary conditions
        NU_R, BB_R, MU_R, WW_R = _flat(m, J)

        # Create subdomains for boundaries
        eps = 1e-6
        axis_test = "x[0] < eps"
        infty_test = "x[0]*x[0] + x[1]*x[1] > (R - eps)*(R - eps)"
        axis  = CompiledSubDomain(axis_test, eps=eps)
        infty = CompiledSubDomain(infty_test, eps=eps, R=R)

        # Create boundary condition on axis for MU
        class AxisValueMU(Expression):
            def eval_cell(self, values, x, cell):
                BB_values = values.copy()
                NU_values = values.copy()
                self.BB.eval_cell(BB_values, x, cell)
                self.NU.eval_cell(NU_values, x, cell)
                values[0] = math.log(BB_values[0]) - NU_values[0]
        MU_a = AxisValueMU(degree=1)
        MU_a.BB = BB
        MU_a.NU = NU

        # Create boundary conditions
        bci_NU = (0, DirichletBC(V, NU_R, infty, method="pointwise"))
        bci_BB = (1, DirichletBC(V, BB_R, infty, method="pointwise"))
        bci_MU = (2, DirichletBC(V, MU_R, infty, method="pointwise"))
        bci_WW = (3, DirichletBC(V, WW_R, infty, method="pointwise"))
        bca_MU = (2, DirichletBC(V, MU_a, axis,  method="pointwise"))
        bc0 = DirichletBC(V, 0.0, DomainBoundary())

        # Collect boundary conditions
        bcs = [bci_NU, bci_BB, bci_MU, bci_WW, bca_MU]

        # Initialize all ansatzes
        for ansatz in ansatzes:
            ansatz.set_fields(NU, BB, MU, WW)
            ansatz.set_integration_parameters(num_steps)
            ansatz.read_parameters()

        # Note: Variables name _foo are unscaled and variables
        # named foo are the corresponding scaled variables:
        # foo = C * _foo.

        # Extract matter fields (unscaled)
        _Phi = model
        _P00 = _Phi[0]
        _P11 = _Phi[1]
        _P33 = _Phi[2]
        _P03 = _Phi[3]

        # Scale matter fields
        C   = Constant(1)
        P00 = C*_P00
        P11 = C*_P11
        P33 = C*_P33
        P03 = C*_P03

        # Define rest density and mass
        rest_density = C*_Phi[4]
        rest_mass    = 2*2*pi*rest_density*r*dx

        # Define density and mass
        _density = BB*(_P00 + _P11 + _P33*(1.0 - (r*BB)**2*WW**2*exp(-4*NU)))
        density  = C*_density
        _mass    = 2*2*pi*_density*r*dx
        mass     = C*_mass

        # Create trial and test functions
        nu = TrialFunction(V)
        bb = TrialFunction(V)
        mu = TrialFunction(V)
        ww = TrialFunction(V)
        v0 = TestFunction(V)
        v1 = TestFunction(V)
        v2 = TestFunction(V)
        v3 = TestFunction(V)

        # Create forms for variational problem
        as_ = _lhs(nu, bb, mu, ww, v0, v1, v2, v3, r, P00, P11, P33, P03)
        Ls  = _rhs(NU, BB, MU, WW, v0, v1, v2, v3, r, P00, P11, P33, P03)

        # Create residuals
        Fs = _lhs(NU, BB, MU, WW, v0, v1, v2, v3, r, P00, P11, P33, P03)
        Fs = [(F - L) for F, L in zip(Fs, Ls)]

        # Create matrices and vectors
        As = [u.vector().factory().create_matrix(mpi_comm_world()) for u in U]
        bs = [u.vector().factory().create_vector(mpi_comm_world()) for u in U]
        Ys = [u.vector().copy() for u in U]

        # Assemble matrices and apply boundary conditions
        [assemble(a, tensor=A) for (a, A) in zip(as_, As)]
        [bc.apply(As[j]) for (j, bc) in bcs]

        # Create linear solver
        preconditioners = [pc for pc in krylov_solver_preconditioners()]
        if "amg" in preconditioners:
            solver = LinearSolver(mpi_comm_world(), "bicgstab", "hypre_amg")            
        else:
            warning("Missing AMG preconditioner, using ILU.")
            solver = LinearSolver(mpi_comm_world(), "bicgstab")

        # Set linear solver parameters
        solver.parameters["relative_tolerance"] = self.parameters.discretization.krylov_tolerance

        # Initialize Anderson acceleration
        anderson = Anderson(depth, [Ui.vector() for Ui in U])

        # Main loop
        tic()
        for iter in xrange(maxiter):

            begin("Iteration %d" % iter)

            # Reset computation of radius of support
            for ansatz in ansatzes:
                ansatz.reset()

            # Rescale right-hand side to preserve mass
            _m = assemble(_mass)
            if _m == 0.0: error("Zero mass distribution.")
            C.assign(m / _m)
            info("C = %.15g" % float(C))

            # Iterate over equations and solve
            for i in range(4):

                # Assemble right-hand side
                assemble(Ls[i], tensor=bs[i])

	        # Apply boundary condition(s)
                [bc.apply(bs[i]) for j, bc in bcs if j == i]

                # Extra left-hand side assembly for B equation
                if i == 1:
                    assemble(as_[i], tensor=As[i])
                    [bc.apply(As[i]) for j, bc in bcs if j == i]

                # Solve linear system
                solver.solve(As[i], Ys[i], bs[i])

                # Damped fixed-point update of solution vector
                #X = U[i].vector()
                #X *= (1.0 - theta)
                #Ys[i] *= theta
                #X += Ys[i]

            # Fixed-point iteration with Anderson acceleration
            Ys = anderson.update_system(Ys)
            for Ui, Yi in zip(U, Ys):
                Ui.vector()[:] = Yi

            # Plot density distribution
            project(density, mesh=mesh, function=RHO)
            self._save_density(RHO, iter)
            self._plot_density(RHO)            

            # Compute residual
            info("Computing residual")
            fs = [assemble(F) for F in Fs]
            [bc0.apply(f) for f in fs]
            residuals = [norm(f, "linf") for f in fs]
            residual = max(r for r in residuals)
            info("||F|| = %.3g (%.3g, %.3g, %.3g, %.3g)" % tuple([residual] + residuals))
            self._save_residual(residual)
            end()

            # Check for convergence
            if residual < tol and iter > 0:
                break

        # Print elapsed time
        info("Iterations finished in %g seconds." % toc())

        # Check whether iteration converged
        solution_converged = False
        if iter == maxiter - 1:
            warning("*** ITERATIONS FAILED TO CONVERGE")
        else:
            info("SOLUTION CONVERGED")
            solution_converged = True
        info("")
        info("Iterations converged to within a tolerance of %g." % tol)
        info("Number of iterations was %g." % iter)
        info("")

        # Compute and store solution characteristics
        self._compute_solution_characteristics(NU, BB, MU, WW, RHO,
                                               P00, P11, P33, P03,
                                               C, _mass, rest_mass,
                                               r, z, V, R,
                                               ansatzes, e0, solution_converged)

        # Compute residuals as functions of space
        fs = [assemble(F) for F in Fs]
        [bc0.apply(f) for f in fs]        
        nu_res = Function(V)
        bb_res = Function(V)
        mu_res = Function(V)
        ww_res = Function(V)
        residual_functions = [nu_res, bb_res, mu_res, ww_res]
        for fres, f in zip(residual_functions, fs):
            fres.vector()[:] = f 

        # Post processing
        solutions = (NU, BB, MU, WW, RHO)
        flat_solutions = (NU_R, BB_R, MU_R, WW_R)
        names = ("NU", "BB", "MU", "WW", "RHO")
        self._postprocess(ansatzes, solutions, flat_solutions, names, residual_functions)

        # Print nice message
        info("Solver complete.")

        return NU, BB, MU, WW, RHO, self.data

    def _compute_solution_characteristics(self,
                                          NU, BB, MU, WW, RHO,
                                          P00, P11, P33, P03,
                                          C, _mass, rest_mass,
                                          r, z, V, R,
                                          ansatzes, e0, solution_converged):
        "Compute interesting properties of solution"

        # Compute final unscaled mass and scale ansatz coefficient
        m = self.parameters.discretization.mass
        _m = assemble(_mass)
        C.assign(m / _m)

        # Compute rest mass
        rm = assemble(rest_mass)

        # Compute total angular momentum
        J = assemble(-2*2*pi*exp(-4.0*NU)*BB*(P03 + (r*BB)**2*WW*P33)*r*dx)

        # Compute gtt component of the metric tensor
        gtt = project(-exp(2.0*NU)*(1.0 - (WW*r*BB)**2*exp(-4.0*NU)), V)

        # Check whether we have an ergo region
        gtt_max = gtt.vector().max()
        ergo_region = gtt_max > 0

        # Compute fractional binding energy
        Eb = 1.0 - m / rm

        # Compute central redshift
        Zc = 1.0/sqrt(-gtt(0,0)) - 1.0

        # Get radius of support and compute areal radius of support
        r0 = max([ansatz.radius_of_support() for ansatz in ansatzes])
        r0 = MPI.max(mpi_comm_world(), r0)
        R0 = r0*(1.0 + m / (2.0*r0))**2

        # Reflection plane radii of support
        rmax = min(2*r0, R)
        rvals = numpy.linspace(0,rmax,10000)
        deltar = rmax/10000.
        RHOvals = numpy.array([RHO(r,0) for r in rvals])
        rp_support = numpy.where(RHOvals > 1e-3)[0] #'vacuum threshold' chosen to agree with 'eye-ball'
        try:
            r_inner = min(rp_support)*deltar
            r_outer = max(rp_support)*deltar
            r_peak  = rvals[numpy.argmax(RHOvals)]
        except ValueError:
            print ("Issue with matter support, setting NA value for rho.")
            r_inner = NA
            r_outer = NA
            r_peak = NA

        # Rcirc
        Rcirc_func = project(r*BB*exp(-NU), V)
        Rcirc = Rcirc_func(r_outer, 0.0)
        
        # Get minimum value of WW
        min_WW = max([ansatz.min_WW() for ansatz in ansatzes])

        # Compute Buchdahl quantity
        Gamma = 2.0*m / R0

        # Store values in dictionary
        self.data = {"ansatz_coefficient": float(C),
                     "mass": m,
                     "radius_of_support": r0,
                     "areal_radius_of_support": R0,
                     "rest_mass": rm,
                     "frac_binding_energy": Eb,
                     "central_redshift": Zc,
                     "gamma": Gamma,
                     "Rcirc": Rcirc,
                     "r_inner": r_inner,
                     "r_outer": r_outer,
                     "r_peak": r_peak,
                     "total_angular_momentum": J,
                     "ergo_region": ergo_region,
                     "gtt_max": gtt_max,
                     "min_WW": min_WW,
                     "solution_converged": solution_converged}

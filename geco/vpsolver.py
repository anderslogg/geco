"""
-----------
vpsolver.py
-----------
This module implements a solver for the Vlasov-Poisson
equations in axial symmetry.

Copyright 2019 Anders Logg, Ellery Ames, Haakan Andreasson

This file is part of GECo.
GECo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GECo is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with GECo. If not, see <https://www.gnu.org/licenses/>.

"""

from solverbase import *


def _lhs(u, v, r):
    return dot(grad(u), grad(v)) * r * dx


def _rhs(rho, v, r):
    return -4 * pi * rho * v * r * dx


def _flat(m):
    return Expression("-m / sqrt(x[0]*x[0] + x[1]*x[1])", degree=1, m=m)


def _init(e0, m, V):
    parameters = {"E": e0, "r0": -m / e0}
    u = "2.0*E / (1.0 + sqrt(x[0]*x[0] + x[1]*x[1]) / r0)"
    return project(Expression(u, degree=1, **parameters), V)


class VlasovPoissonSolver(SolverBase):
    "Solver for the axisymmetric Vlasov-Poisson equations"

    def __init__(self):
        SolverBase.__init__(self, "vp")

    def solve(self, model, solution=None):
        "Compute solution"

        # Initialize solve
        self._init_solve()

        # Extract all ansatzes from model
        ansatzes = [
            c for c in extract_coefficients(model) if not isinstance(c, Constant)
        ]

        # Get common model parameters (use first)
        e0 = ansatzes[0].parameters.E0

        # Get discretization parameters
        m = self.parameters.discretization.mass
        R = self.parameters.discretization.domain_radius
        maxiter = self.parameters.discretization.maxiter
        tol = self.parameters.discretization.tolerance
        num_steps = self.parameters.discretization.num_steps
        degree = self.parameters.discretization.degree
        depth = self.parameters.discretization.anderson_depth

        # Get output parameters
        plot_iteration = self.parameters.output.plot_iteration

        # Workaround for geometric round-off errors
        parameters.allow_extrapolation = True

        # Generate mesh and create function space
        mesh = self._generate_mesh()
        V = FunctionSpace(mesh, "Lagrange", degree)

        # Create initial data
        U = _init(e0, m, V)
        RHO = Function(V)

        # Override initial data if provided
        if solution is not None:
            info("Reusing function space and initial data.")
            U, RHO = solution[:2]
            V = U.function_space()
            mesh = V.mesh()

        # Create spatial coordinates
        x = SpatialCoordinate(mesh)
        r = x[0]
        z = x[1]

        # Create asymptotically flat solutions for boundary conditions
        U_R = _flat(m)

        # Create subdomains for boundaries
        eps = 1e-3
        infty_test = "x[0]*x[0] + x[1]*x[1] > (R - eps)*(R - eps)"
        infty = CompiledSubDomain(infty_test, eps=eps, R=R)

        # Create boundary conditions
        bc = DirichletBC(V, U_R, infty)
        bc0 = DirichletBC(V, 0.0, infty)

        # Initialize all ansatzes
        for ansatz in ansatzes:
            ansatz.set_fields(U)
            ansatz.set_integration_parameters(num_steps)
            ansatz.read_parameters()

        # Note: Variables name _foo are unscaled and variables
        # named foo are the corresponding scaled variables:
        # foo = C * _foo.

        # Extract density (unscaled)
        _density = model

        # Scale density
        C = Constant(1)
        density = C * _density

        # Define unscaled and scaled mass
        _mass = 2 * 2 * pi * _density * r * dx
        mass = C * _mass

        # Create trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)

        # Create forms for variational problem
        a = _lhs(u, v, r)
        L = _rhs(density, v, r)

        # Create residual
        F = _lhs(U, v, r) - L

        # Assemble matrix and apply boundary condition
        A = assemble(a)
        bc.apply(A)

        # Create vector
        b = Vector()
        Y = U.vector().copy()

        # Create linear solver
        preconditioners = [pc for pc in krylov_solver_preconditioners()]
        if "amg" in preconditioners:
            solver = LinearSolver(mpi_comm_world(), "gmres", "amg")
        else:
            warning("Missing AMG preconditioner, using ILU.")
            solver = LinearSolver(mpi_comm_world(), "gmres")

        # Set linear solver parameters
        solver.parameters[
            "relative_tolerance"
        ] = self.parameters.discretization.krylov_tolerance

        # Initialize Anderson acceleration
        anderson = Anderson(depth, U.vector())

        # Main loop
        tic()
        for iter in xrange(maxiter):

            begin("Iteration %d" % iter)

            # Reset computation of radius of support
            for ansatz in ansatzes:
                ansatz.reset()

            # Rescale right-hand side to preserve mass
            _m = assemble(_mass)
            if _m == 0.0:
                error("Zero mass distribution.")
            C.assign(m / _m)
            info("C = %.15g" % float(C))

            # Assemble right-hand side and apply boundary conditions
            assemble(L, tensor=b)
            bc.apply(b)

            # Solve linear system
            solver.solve(A, Y, b)

            # Fixed-point iteration with Anderson acceleration
            U.vector()[:] = anderson.update(Y)

            # Damped fixed-point update of solution vector
            # theta = self._get_theta()
            # X = U.vector()
            # X *= (1.0 - theta)
            # Y *= theta
            # X += Y

            # Plot density distribution
            project(density, mesh=mesh, function=RHO)
            self._save_density(RHO, iter)
            self._plot_density(RHO)

            # Compute residual
            info("Computing residual")
            f = assemble(F)
            bc0.apply(f)
            residual = norm(f, "linf")
            info("||F|| = %.3g" % residual)
            self._save_residual(residual)
            end()

            # Check for convergence
            if residual < tol and iter > 0:
                break

        # Print elapsed time
        info("Iterations finished in %g seconds." % toc())

        # Check whether iteration converged
        if iter == maxiter - 1:
            warning("*** ITERATIONS FAILED TO CONVERGE")
        else:
            info("SOLUTION CONVERGED")
        info("")
        info("Iterations converged to within a tolerance of %g." % tol)
        info("Number of iterations was %g." % iter)
        info("")

        # Compute and store solution characteristics
        self._compute_solution_characteristics(C, _mass, ansatzes)

        # Compute residuals as functions of space
        fs = assemble(F)
        bc0.apply(fs)
        U_res = Function(V)
        U_res.vector()[:] = fs
        residual_functions = (U_res,)

        # Post processing
        solutions = (U,)
        flat_solutions = (U_R,)
        matter_components = (RHO,)
        names = ("U",)
        matter_names = ("RHO",)
        # self._postprocess(ansatzes, solutions, flat_solutions, names)
        self._postprocess(
            ansatzes,
            solutions,
            flat_solutions,
            names,
            residual_functions,
            matter_components,
            matter_names,
        )

        # Print nice message
        info("Solver complete")

        return U, RHO, self.data

    def _compute_solution_characteristics(self, C, _mass, ansatzes):
        "Compute interestng properties of solution"

        # Compute final unscaled mass and scale ansatz coefficient
        m = self.parameters.discretization.mass
        _m = assemble(_mass)
        C.assign(m / _m)

        # Get radius of support and compute areal radius of support
        r0 = max([ansatz.radius_of_support() for ansatz in ansatzes])
        r0 = MPI.max(mpi_comm_world(), r0)
        R0 = r0 * (1.0 + m / (2.0 * r0)) ** 2

        # Compute Buchdahl quantity
        Gamma = 2.0 * m / R0

        # Store results
        self.data = {
            "ansatz_coefficient": float(C),
            "mass": m,
            "radius_of_support": r0,
            "areal_radius_of_support": R0,
            "gamma": Gamma,
        }

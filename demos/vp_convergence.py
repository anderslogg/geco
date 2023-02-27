"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Polytropic

This creates a torus solution. This demo is similar to the vp_torus
demo but additionally sets a much stricter tolerance and monitors the
convergence of the iterations. Use geco-postprocess-convergence to
plot the residual convergence.
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters.discretization.tolerance = 1e-10
solver.parameters.discretization.krylov_tolerance = 1e-16
solver.parameters.discretization.anderson_depth = 9

# Create ansatz
model = MaterialModel("VP-E-Polytropic-L-Polytropic")
model.parameters.L0 = 1.0
model.parameters.l = 1.0

# Solve
U, RHO, data = solver.solve(model)

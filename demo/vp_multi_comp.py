"""
Equation: Vlasov-Poisson
Ansatz:   multi-species test
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.domain_radius = 50
solver.parameters.discretization.resolution = 32
solver.parameters.discretization.tolerance = 1e-3

# Create ansatzes
S = MaterialModel("VP-E-Polytropic-L-Gaussian")
T = MaterialModel("VP-E-Polytropic-L-Polytropic")
U = MaterialModel("VP-E-Polytropic-L-Polytropic")
model = Constant(0.3)*T + Constant(0.3)*S + Constant(0.4)*U

# Ansatz parameters
S.parameters.E0 = -0.06
S.parameters.L0 = 1.1
S.parameters.k = 2.4
S.parameters.pm = 1.0

T.parameters.E0 = -0.1
T.parameters.L0 = 1.0
T.parameters.k = 2.0

U.parameters.E0 = -0.2
U.parameters.L0 = 1.0
U.parameters.k = 0.5

# Compute solution for initial guess
solution = solver.solve(model)

# Extract solution components
#NU, BB, MU, WW, RHO, data = solution

# Plot density
#plot(RHO)
#interactive()

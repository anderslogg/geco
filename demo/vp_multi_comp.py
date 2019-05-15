"""
Equation: Vlasov-Poisson
Ansatz:   multi-species test
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.domain_radius = 50
solver.parameters.discretization.resolution = 256
solver.parameters.discretization.tolerance = 1e-6

# Create ansatzes
S = MaterialModel("VP-E-Polytropic-L-Gaussian")
T = MaterialModel("VP-E-Polytropic-L-Polytropic")
model = S #Constant(0.5)*T #Constant(0.5)*S + 

# Ansatz parameters
S.parameters.E0 = -0.06
S.parameters.L0 = 1.1
S.parameters.k = 2.4
S.parameters.pm = 1.0

# Compute solution for initial guess
solution = solver.solve(model)

# Extract solution components
#NU, BB, MU, WW, RHO, data = solution

# Plot density
#plot(RHO)
#interactive()

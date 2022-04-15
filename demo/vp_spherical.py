"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Polytropic

This creates a spherical solution of the VP system. 

This demo converges in two iterations. 
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters.discretization.domain_radius = 30
solver.parameters.discretization.resolution = 256
solver.parameters.discretization.tolerance = 1e-8

# Create ansatz
model = MaterialModel("VP-E-Polytropic-L-Polytropic")
model.parameters.L0 = 0.0
model.parameters.E0 = -0.05
model.parameters.l = 0.0

# Solve
U, RHO, data = solver.solve(model)

"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Polytropic

This creates a spherical solution of the VP system. 

This demo converges in four iterations. 
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters["discretization"]["domain_radius"] = 30
solver.parameters["discretization"]["resolution"]    = 64
solver.parameters["discretization"]["tolerance"]     = 1e-6

# Create ansatz
model = MaterialModel("VP-E-Polytropic-L-Andreasson")
model.parameters["E0"] = -0.1
model.parameters["k"]  = 0.0
model.parameters["Q"] = 0.1
model.parameters["l"]  = 0.0

# Solve
U, RHO, data = solver.solve(model)

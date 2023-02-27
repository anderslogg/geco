"""
Equation: Vlasov-Poisson
Ansatz:   VP-Rowley

This demo converges in four iterations. 
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters["discretization"]["domain_radius"] = 30
solver.parameters["discretization"]["resolution"]    = 64
solver.parameters["discretization"]["tolerance"]     = 1e-6

# Create ansatz
model = MaterialModel("VP-Rowley")
model.parameters["E0"] = -0.1
model.parameters["r0"] = 1.0
model.parameters["W"]  = 1.0
model.parameters["c0"] = -0.1
model.parameters["b"]  = 0.1

# Solve
U, RHO, data = solver.solve(model)

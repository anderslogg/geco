"""
Equation: Vlasov-Poisson
Ansatz:   VP-Evans-L-Polytropic

This demo converges in seven iterations. 
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters["discretization"]["domain_radius"] = 30
solver.parameters["discretization"]["resolution"]    = 64
solver.parameters["discretization"]["tolerance"]     = 1e-6

# Create ansatz
model = MaterialModel("VP-Evans-L-Polytropic")
model.parameters["E0"] = -0.1
model.parameters["v0"]  = 0.0
model.parameters["s0"] = 1.0
model.parameters["c"]  = 0.1

# Solve
U, RHO, data = solver.solve(model)

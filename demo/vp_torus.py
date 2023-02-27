"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Polytropic

This creates a solution of the VP system with toroidal morphology.

Converges in ten iterations.
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters["discretization"]["tolerance"]= 1e-6
solver.parameters["discretization"]["resolution"] = 64

# Create ansatz
model = MaterialModel("VP-E-Polytropic-L-Polytropic")
model.parameters["L0"] = 1.0
model.parameters["l"]= 1.0

# Solve
U, RHO, data = solver.solve(model)

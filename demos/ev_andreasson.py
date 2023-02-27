"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Andreasson

This creates a non-rotating EV solution with spindle-like morphology.

The final solution converges in 3 iterations.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters["output"]["plot_solution"] = False
solver.parameters["discretization"]["domain_radius"] = 30
solver.parameters["discretization"]["resolution"] = 64

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Andreasson")
model.parameters["rotation"] = False

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model.parameters["Q"] = 2.5

# Compute solution
solution = solver.solve(model, solution)

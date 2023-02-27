"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Gaussian

This creates an EV solution with prolate spheroid (spindle) morphology and zero net angular momentum.

Final solution converges in 2 iterations.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters["output"]["plot_solution"] = False
solver.parameters["discretization"]["domain_radius"] = 100
solver.parameters["discretization"]["resolution"] = 128

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters["E0"] = 0.966
model.parameters["k"]= 0.0
model.parameters["l"]= 0.0
model.parameters["rotation"] = False

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("EV-E-Polytropic-L-Gaussian")
model.parameters["E0"] = 0.966
model.parameters["pm"] = -1.0
model.parameters["k"] = 0.0
model.parameters["L0"] = 0.1
model.parameters["rotation"] = False

# Compute solution
solution = solver.solve(model, solution)

"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Gaussian

This creates a self-gravitating solution with oblate spheroidal morphology.
Two ansatzes are used, one to generate a self-gravitating ball, which is then fed as an initial guess for the solver with the disk ansatz. 
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters["output"]["plot_solution"] = False
solver.parameters["discretization"]["domain_radius"] = 50

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters["E0"] = 0.942
model.parameters["k"]= 0.0
model.parameters["l"]= 0.0
model.parameters["rotation"] = False

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("EV-E-Polytropic-L-Gaussian")
model.parameters["E0"] = 0.942
model.parameters["pm"] = 1.0
model.parameters["k"] = 2.0
model.parameters["L0"] = 1.4
model.parameters["rotation"] = False

# Compute solution
solution = solver.solve(model, solution)

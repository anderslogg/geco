"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Gaussian

This creates a rotating solution of the EV system in which all particles have nearly identical energy and angular momentum. 
The morphology is that of a torus, with energy density concentrated at the inner edge. 
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters["output"]["plot_solution"] = False
solver.parameters["discretization"]["domain_radius"] = 50
solver.parameters["discretization"]["resolution"] = 128
solver.parameters["discretization"]["num_steps"] = 32

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters["E0"] = 0.92
model.parameters["k"] = 0.0
model.parameters["l"] = 0.0
model.parameters["particle_mass"] = 1.0

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("EV-E-Delta-L-Delta")
model.parameters["E0"] = 0.92
model.parameters["L0"] = 0.5

# Compute solution
solution = solver.solve(model, solution)

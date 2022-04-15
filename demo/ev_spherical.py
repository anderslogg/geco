"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates a non-rotating spherically symmetric solution.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.925
model.parameters.rotation = False

# Compute solution 
solution = solver.solve(model)

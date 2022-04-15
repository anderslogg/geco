"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates a rotating solution of the EV system with toroidal morphology.

To converge, this demo should be run at higher resolution. 
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False

# If not in hires, increase resolution
solver.parameters.discretization.num_steps = 32
solver.parameters.discretization.resolution = 64

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz (adjust L0 from 0 to 1)
model.parameters.L0 = 1.0

# Compute solution
solution = solver.solve(model, solution)

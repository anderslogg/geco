"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates a torus solution.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.radius = 50

#solver.parameters.discretization.resolution = 64

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.925
model.parameters.rotation = False

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz (adjust L0 from 0 to 1)
model.parameters.k = 1.0
model.parameters.l = 1.0

# Compute solution
solution = solver.solve(model, solution)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Plot density
plot(RHO)
interactive()

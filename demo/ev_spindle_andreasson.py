"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Andreasson

This creates a non-rotating EV spindle solution.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.radius = 50

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Andreasson")
model.parameters.rotation = False

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model.parameters.Q = 2.5

# Compute solution
solution = solver.solve(model, solution)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Plot density
plot(RHO)
interactive()

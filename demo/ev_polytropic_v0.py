"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates a non-rotating torus solution.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.radius = 50

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.925
model.parameters.rotation = False

# Compute solution 
solution = solver.solve(model)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Plot density
plot(RHO)
interactive()

"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Gaussian

This creates an EV spindle solution with zero net angular momentum.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.radius = 100

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.966
model.parameters.k = 0.0
model.parameters.l = 0.0
model.parameters.rotation = False

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("EV-E-Polytropic-L-Gaussian")
model.parameters.E0 = 0.966
model.parameters.pm = -1.0
model.parameters.k  = 0.0
model.parameters.L0 = 0.1
model.parameters.rotation = False

# Compute solution
solution = solver.solve(model, solution)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Plot density
plot(RHO)
interactive()

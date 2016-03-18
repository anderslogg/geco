"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Gaussian

This creates a rotating disk of particles with positive angular momentum.

Note: Because of the large E0 parameter, this dem should be run on a mesh of at least 50.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.radius = 50
solver.parameters.discretization.resolution = 64

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.942
model.parameters.k = 0.0
model.parameters.l = 0.0

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("EV-E-Polytropic-L-Gaussian")
model.parameters.E0 = 0.942
model.parameters.pm = 1.0
model.parameters.k  = 1.6
model.parameters.L0 = 1.27

# Compute solution
solution = solver.solve(model, solution)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Plot density
plot(RHO)
interactive()

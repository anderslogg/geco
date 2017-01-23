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
solver.parameters.discretization.num_steps = 32

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.92
model.parameters.k = 0.0
model.parameters.l = 0.0
model.parameters.particle_mass = 1.0

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("EV-E-Delta-L-Delta")
model.parameters.E0 = 0.92
model.parameters.L0 = 0.5

# Compute solution
solution = solver.solve(model, solution)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Check results
J = data["total_angular_momentum"]
M = data["rest_mass"]
m = data["particle_mass"]
eps = abs(J - M/m*0.5)
print eps

"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates a torus solution.

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

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Check results
J = data["total_angular_momentum"]
M = data["rest_mass"]
m = data["particle_mass"]
eps = abs(J - M/m)
print eps

# Plot density
plot(RHO)
interactive()

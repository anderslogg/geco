"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates a thin relativistic torus solution with a low E0 value.
In comparison to the EVR version this demo demonstrates that rotation stabilizes
the solutions.

"""

from geco import *
from numpy import linspace
from numpy import floor
import sys

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False

# Set parameters
solver.parameters.discretization.tolerance = 1e-4
solver.parameters.discretization.radius = 25
#solver.parameters.discretization.num_steps = 25
solver.parameters.discretization.theta = 1.0

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.rotation = False

# Compute solution for initial guess
solution = solver.solve(model)

# Adjust value of L0
model.parameters.L0 = 0.5
solution = solver.solve(model, solution)

# Step down the value of E0
e0vals = list(linspace(0.90, 0.82, 10))
for e in e0vals:
    print "E0 =", e
    model.parameters.E0 = e
    solution = solver.solve(model, solution)

# Decrease tolerance and theta to step a little more carefully
solver.parameters.discretization.tolerance = 1e-4
solver.parameters.discretization.theta = 0.5

# Step down the value of E0
e0vals = linspace(0.82, 0.78, 10)
for e in e0vals:
    print "E0 =", e
    model.parameters.E0 = e
    solution = solver.solve(model, solution)

# Decrease theta further
solver.parameters.discretization.theta = 0.25

# Step down the value of E0
e0vals = linspace(0.78, 0.70, 20)
for e in e0vals:
    print "E0 =", e
    model.parameters.E0 = e
    solution = solver.solve(model, solution)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Plot density
plot(RHO)
interactive()

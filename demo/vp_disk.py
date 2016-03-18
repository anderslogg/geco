"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Gaussian

This creates a "thin" extremal disk solution. By extremal we mean
that increasing k further makes the body too compact and un-disk
like, while decreasing L0 further produces a thin torus at large
radius.
"""

from geco import *
from numpy import linspace

# Create solver
solver = VlasovPoissonSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.radius = 50

# Create ansatz for initial guess
model = MaterialModel("VP-E-Polytropic-L-Polytropic")
model.parameters.E0 = -0.06
model.parameters.L0 = 0.0
model.parameters.k = 1.0
model.parameters.l = 1.0

# Compute solution for initial guess
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("VP-E-Polytropic-L-Gaussian")
model.parameters.E0 = -0.06
model.parameters.pm = 1.0
model.parameters.k = 2.4

# Gradually decrease L0
for L0 in list(linspace(1.4, 1.2, 3)) + list(linspace(1.2, 1.1, 3)):
    print "L0 =", L0

    # Adjust L0
    model.parameters.L0 = L0

    # Compute solution
    solution = solver.solve(model, solution)

# Extract solution components
U, RHO, data = solution

# Plot density
plot(RHO)
interactive()

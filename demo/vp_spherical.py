"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Polytropic

This creates a spherical solution.
"""

from geco import *

# Create solver
solver = VlasovPoissonSolver()
solver.parameters.output.plot_solution = False

# Create ansatz
model = MaterialModel("VP-E-Polytropic-L-Polytropic")
model.parameters.L0 = 0.0
model.parameters.l = 0.0

# Solve
U, RHO, data = solver.solve(model)

# Plot density
plot(RHO)
interactive()

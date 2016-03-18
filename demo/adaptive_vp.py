"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Polytropic

This is a first simple test of the adaptivity. Nothing special
or fancy about the solution itself.
"""

from geco import *

from fenics import *
parameters.reorder_dofs_serial = False

# Create solver
solver = AdaptiveVlasovPoissonSolver()

# Create ansatz
model = MaterialModel("VP-E-Polytropic-L-Polytropic")

# Compute solution
U, RHO, data = solver.solve(model)

# Plot density
#plot(RHO)
#interactive()

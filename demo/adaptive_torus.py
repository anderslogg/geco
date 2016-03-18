"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Polytropic

This is a simple test of the adaptivity to create VP torus.

Note the code converges for L0 = 0.1, but not greater. 
"""

from geco import *

from fenics import *
parameters.reorder_dofs_serial = False

# Create solver
solver = AdaptiveVlasovPoissonSolver()

# Create ansatz
model = MaterialModel("VP-E-Polytropic-L-Polytropic")
model.parameters.L0 = 1.0

# Compute solution
U, RHO, data = solver.solve(model)

# Plot density
#plot(RHO)
#interactive()

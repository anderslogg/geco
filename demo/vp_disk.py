"""
Equation: Vlasov-Poisson
Ansatz:   VP-E-Polytropic-L-Gaussian

This creates solution of the VP system with a oblate (disk) morphology.
The demo also demonstrates how different solutions can be saved in different files, and how to pass a prior solution as an initial guess for a subsequent ansatz.
"""

from geco import *
from numpy import linspace

# Output directory
out_dir = "vp_oblate_solutions/"

# Create solver
solver = VlasovPoissonSolver()
solver.parameters["discretization"]["domain_radius"] = 50

# Create ansatz for initial guess
model = MaterialModel("VP-E-Polytropic-L-Polytropic")
model.parameters["E0"] = -0.06
#model.parameters.E0 = -0.06
model.parameters.L0 = 0.0
model.parameters.k = 1.0
model.parameters.l = 1.0

# Compute solution for initial guess
solver.parameters.output.solution_directory = out_dir + 'init_ball'
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("VP-E-Polytropic-L-Gaussian")
model.parameters.E0 = -0.06
model.parameters.pm = 1.0
model.parameters.k = 2.4

# Gradually decrease L0
for L0 in list(linspace(1.4, 1.2, 3)) + list(linspace(1.2, 1.1, 3)):
    print('L0 =', L0)

    # Change solution output directory
    solver.parameters.output.solution_directory = out_dir + 'L0_{}'.format(L0)

    # Adjust L0
    model.parameters.L0 = L0

    # Compute solution
    solution = solver.solve(model, solution)

# Extract solution components
U, RHO, data = solution

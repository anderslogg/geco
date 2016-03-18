"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Andreasson (spindle)
          EV-E-Polytropic-L-Polytropic (torus) 

This creates an L0-parametrized family of hoag-type objects with a spindle center and a toroidal ring.

Each has zero net angular momentum.

For L0 larger than 1.8, the density lies almost all in the spindle component.

"""

from geco import *
import os
from numpy import linspace

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
#solver.parameters.discretization.radius = 200
solution_dir = solver.parameters.output.solution_directory

# Create ansatzes
S = MaterialModel("EV-E-Polytropic-L-Andreasson")
T = MaterialModel("EV-E-Polytropic-L-Polytropic")
model = Constant(0.5)*S + Constant(1.0)*T

# Ansatz parameters
S.parameters.E0 = 0.940
S.parameters.Q  = 2.0
S.parameters.k  = 1.0
S.parameters.l  = 0.0
S.parameters.rotation = False

T.parameters.E0 = 0.940
T.parameters.L0 = 0.0
T.parameters.k = 0.5
T.parameters.l = 1.0
T.parameters.rotation = False

# Family of L0-parametrized solutions
L0vals = list(linspace(1.3, 1.8, 11))
for l in L0vals:
    print "L0 =", l
    T.parameters.L0 = l

    # Set up directories for solution and data
    specificL_dir = os.path.join(solution_dir,"L0_%g" % l)
    solver.parameters.output.solution_directory = specificL_dir

    # Solve (starting from default initial guess)
    solution = solver.solve(model)


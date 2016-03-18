"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Polytropic

The total mass is varied, as is the parameter L0, in order to verify scaling of solutions.

"""

from geco import *
from numpy import linspace
from numpy import floor
import sys, os

# The change_dir function
def change_dir(m):
    # Set name / prefix of simulation
    specificM_dir = os.path.join(solution_dir,"Mass_%g" % m)
    solver.parameters.output.solution_directory = specificM_dir
    
# Create solver
solver = EinsteinVlasovSolver()
solution_dir = solver.parameters.output.solution_directory

# Set parameters
prefix = "halfcircle"
solver.parameters.discretization.tolerance = 1e-4
solver.parameters.discretization.radius = 100
solver.parameters.discretization.resolution = 32
solver.parameters.discretization.num_steps = 32
solver.parameters.discretization.theta = 1.0
solver.parameters.discretization.mesh_prefix = prefix
solver.parameters.output.plot_solution = False

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")

# Compute base solution
solver.parameters.discretization.mass = 1.0
model.parameters.L0 = 0.0
model.parameters.k  = 0.0
model.parameters.l  = 0.0
m = 1.0
change_dir(m)
solution = solver.solve(model)

# Vary the total mass
massvals = list(linspace(1.0, 5.0, 17))
for m in massvals:
    print "Mass =", m
    change_dir(m)
    solver.parameters.discretization.mass = m
    solution = solver.solve(model, solution)


"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Gaussian

This creates a family of flattened spheroids. 

This demo uses a non-uniform mesh and should NOT be run in hires mode. 
"""

from geco import *
from numpy import linspace
from numpy import floor
import sys, os

# The change_dir function
def change_dir(l):
    # Set name / prefix of simulation
    specificL_dir = os.path.join(solution_dir,"L0_%g" % l)
    solver.parameters.output.solution_directory = specificL_dir
    
# Create solver
solver = EinsteinVlasovSolver()

# Configure solution directory if not provided
solution_dir = solver.parameters.output.solution_directory   

# Set parameters
prefix = "nested_halfdisk"
solver.parameters.discretization.tolerance = 1e-2
solver.parameters.discretization.radius = 100
solver.parameters.discretization.resolution = 16
solver.parameters.discretization.num_steps = 25
solver.parameters.discretization.theta = 1.0
solver.parameters.discretization.mesh_prefix = prefix
solver.parameters.output.plot_solution = False

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.942
model.parameters.k = 0.0
model.parameters.l = 0.0
model.parameters.rotation = False

# Compute solution for initial guess
L0 = 10.0
change_dir(L0)
solution = solver.solve(model)

# Create main ansatz
model = MaterialModel("EV-E-Polytropic-L-Gaussian")
model.parameters.E0 = 0.942
model.parameters.pm = 1.0
model.parameters.k  = 1.5
model.parameters.L0 = 10.0
model.parameters.rotation = False

# Compute solution for main ansatz
solver.parameters.discretization.tolerance = 1e-4
solution = solver.solve(model, solution)

# Step up the value of L0
l0vals = list(linspace(10.0, 3.0, 5)) + list(linspace(3.0, 2.1, 10))
for l in l0vals:
    print "L0 =", l
    model.parameters.L0 = l
    # If L0 in saved_L0, change solution directory
    if l not in l0vals[0:16:4]:
        solution = solver.solve(model, solution)
    else:
        change_dir(l)
        solution = solver.solve(model, solution)

# Step up the value of L0
l0vals = linspace(2.0, 1.0, 11)
for l in l0vals:
    print "L0 =", l
    model.parameters.L0 = l
    change_dir(l)
    solution = solver.solve(model, solution)        

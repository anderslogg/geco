"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates a thin relativistic torus solution.
The angular momentum is less than for the L0 = 0.8 version
and eventually J/M^2 < 1

This demo uses a non-uniform mesh and should NOT be run in hires mode.
mesh: ccr_halfdisk_50_16

"""

from geco import *
from numpy import linspace
from numpy import floor
import sys, os

# L0 Parameter
L0 = 0.5

# The change_dir function
def change_dir(e):
    # Set name / prefix of simulation
    specificE_dir = os.path.join(solution_dir,"E0_%g" % e)
    solver.parameters.output.solution_directory = specificE_dir
    
# Create solver
solver = EinsteinVlasovSolver()

# Configure solution directory if not provided
solution_dir = solver.parameters.output.solution_directory
if solution_dir == "solutions/evr":
    solution_dir = os.path.join("solutions", "evr_rel_torus_%g" % L0)
    

# Set parameters
prefix = "ccr_halfdisk"
solver.parameters.discretization.tolerance = 1e-3
solver.parameters.discretization.radius = 50
solver.parameters.discretization.resolution = 16
solver.parameters.discretization.num_steps = 32
solver.parameters.discretization.theta = 0.75
solver.parameters.discretization.mesh_prefix = prefix
solver.parameters.output.plot_solution = False

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")

# Compute solution for initial guess
e = 0.9
change_dir(e)
solution = solver.solve(model)

# Adjust value of L0
model.parameters.L0 = L0
solution = solver.solve(model, solution)

# Step down the value of E0
e0vals = list(linspace(0.90, 0.80, 10)) + list(linspace(0.80, 0.77, 10))
for e in e0vals:
    print "E0 =", e
    model.parameters.E0 = e
    # If E0 in saved_E0, change solution directory
    if e not in e0vals[0:20:5]:
        solution = solver.solve(model, solution)
    else:
        change_dir(e)
        solution = solver.solve(model, solution)

# Decrease tolerance and theta to step a little more carefully
solver.parameters.discretization.tolerance = 1e-3
solver.parameters.discretization.theta = 0.5

# Step down the value of E0
e0vals = linspace(0.77, 0.75, 10)
for e in e0vals:
    print "E0 =", e
    model.parameters.E0 = e
    # If E0 in saved_E0, change solution directory
    if e not in e0vals[0:10:5]:
        solution = solver.solve(model, solution)
    else:
        change_dir(e)
        solution = solver.solve(model, solution)

# Decrease theta further
solver.parameters.discretization.tolerance = 1e-4
solver.parameters.discretization.theta = 0.25

# Step down the value of E0
e0vals = linspace(0.75, 0.55, 101)
for e in e0vals:
    print "E0 =", e
    model.parameters.E0 = e
    # If E0 in saved_E0, change solution directory
    if e not in e0vals[0:101:5]:
        solution = solver.solve(model, solution)
    else:
        change_dir(e)
        solution = solver.solve(model, solution)

# Keep stepping down the value of E0 if possible
e0vals = linspace(0.55, 0.45, 51)
for e in e0vals:
    print "E0 =", e
    model.parameters.E0 = e
    # If E0 in saved_E0, change solution directory
    if e not in e0vals[0:51:5]:
        solution = solver.solve(model, solution)
    else:
        change_dir(e)
        solution = solver.solve(model, solution)
        

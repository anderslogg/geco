"""
Equation: Vlasov-Poisson, Einstein-Vlasov, Einstein-Vlasov-Rotating
Ansatz:   E-Polytropic-L-Gaussian

This creates a disk-like solution with the same
parameters in each model mentioned above.

This demo uses a non-uniform mesh and should NOT be run in hires mode. 
"""

from geco import *
from numpy import linspace
from numpy import floor
import sys, os


solvers  = [VlasovPoissonSolver(), EinsteinVlasovSolver(), EinsteinVlasovSolver()]
ansatzes = ["VP-E-Polytropic-L-Gaussian", "EV-E-Polytropic-L-Gaussian","EV-E-Polytropic-L-Gaussian"]
E0params = [-0.06, 0.942, 0.942]
names    = ["vp","ev","evr"]

def _setup_model(solver, ansatz, e0, name):

    # Set parameters
    prefix = "nested_halfdisk"
    solver.parameters.discretization.tolerance = 1e-2
    solver.parameters.discretization.radius = 100
    solver.parameters.discretization.resolution = 16
    solver.parameters.discretization.num_steps = 32
    solver.parameters.discretization.theta = 1.0
    solver.parameters.discretization.mesh_prefix = prefix
    solver.parameters.output.plot_solution = False

    # Set solution directory
    solution_dir = solver.parameters.output.solution_directory  
    solver.parameters.output.solution_directory = os.path.join(solution_dir,name + "_disk")

    # Create ansatz for initial guess
    model = MaterialModel(ansatz)
    model.parameters.E0 = e0
    model.parameters.k  = 0.0
    model.parameters.L0 = 10.0
    model.parameters.pm = 1.0

    # Handle rotation in EV
    if name == "ev":
        model.parameters.rotation = False

    # Return model
    return model    


for (solver,ansatz,e0,name) in zip(solvers,ansatzes,E0params,names):

    # Print message
    info("Computing disk solution in %s model" % name)

    # Set parameters
    model = _setup_model(solver, ansatz, e0, name)

    # Compute initial guess
    solution = solver.solve(model)

    # Set main ansatz parameters
    model.parameters.k  = 1.5
    model.parameters.L0 = 1.5

    # Compute solution
    solver.parameters.discretization.tolerance = 1e-4
    solution = solver.solve(model,solution)
    


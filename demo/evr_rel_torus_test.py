"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Polytropic

This is a TEST evr rel torus demo to
understand the effects of damping on
the ripples in the numerical solution.

This demo uses a non-uniform mesh and should NOT be run in hires mode.

"""

from geco import *
from numpy import linspace
from numpy import floor
import sys, os

# Use previous initial data if any
use_initial_data = True

# Create solver
solver = EinsteinVlasovSolver()

# Set parameters
prefix = "nested_halfdisk"
solver.parameters.discretization.tolerance = 1e-3
solver.parameters.discretization.radius = 50
solver.parameters.discretization.resolution = 32
solver.parameters.discretization.num_steps = 25
solver.parameters.discretization.theta = 1.0
solver.parameters.discretization.mesh_prefix = prefix
solver.parameters.output.plot_solution = False

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")

# Either iterate or start from precomputed initial state
if use_initial_data:

    # Create mesh and initialize solutions
    mesh = Mesh("../geco/geometries/nested_halfdisk_50_32_mesh.xml.gz")
    V = FunctionSpace(mesh, "Lagrange", 1)
    NU = Function(V)
    BB = Function(V)
    MU = Function(V)
    WW = Function(V)
    RHO = Function(V)
    File("evr_rel_torus_test_data/NU_50.xml") >> NU
    File("evr_rel_torus_test_data/BB_50.xml") >> BB
    File("evr_rel_torus_test_data/MU_50.xml") >> MU
    File("evr_rel_torus_test_data/WW_50.xml") >> WW
    File("evr_rel_torus_test_data/RHO_50.xml") >> RHO

    # Compute solution
    model.parameters.L0 = 0.8
    model.parameters.E0 = 0.8
    solution = solver.solve(model, (NU, BB, MU, WW, RHO))

else:

    # Compute solution for initial guess
    solution = solver.solve(model)

    # Adjust value of L0
    model.parameters.L0 = 0.8
    solution = solver.solve(model, solution)

    # Step down the value of E0
    e0vals = list(linspace(0.90, 0.80, 10))
    for e in e0vals:
        print "E0 =", e
        model.parameters.E0 = e
        solution = solver.solve(model, solution)

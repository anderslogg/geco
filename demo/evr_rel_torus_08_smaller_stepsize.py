"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Polytropic

This creates our most relativistic thin torus solution.

This demo recomputes the solution starting from the final iteration
computed with evr_rel_torus_08.py with a smaller step size for the
numerical integration to get rid of the ripples in the solution.

Relevant values for E0:

  0.58  0.6   0.62  0.64  0.66  0.68  0.7   0.72
  0.59  0.61  0.63  0.65  0.67  0.69  0.71  0.73

"""

from geco import *
import sys

# Get value of E0
E0 = 0.58
## if not len(sys.argv) == 2:
##     error("Missing argument E0.")
## E0 = sys.argv[1]

# Set data directory
data_dir = "/home/ellery/solutions/paper/evr_rel_torus_midL-2016-03-04-17-10-50/E0_"
data_dir = data_dir + "0.58"

# Create model
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.L0 = 0.8

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.discretization.radius = 50
solver.parameters.discretization.resolution = 16
solver.parameters.discretization.num_steps = 100
solver.parameters.discretization.mesh_prefix = "ccr_halfdisk"
solver.parameters.discretization.tolerance = 1e-4
solver.parameters.discretization.theta = 0.25
solver.parameters.output.plot_solution = False
solver.parameters.output.solution_directory = "solutions/" + "{:.9f}".format(E0) 

# Read initial data
mesh = Mesh("../geco/geometries/ccr_halfdisk_50_16_mesh.xml.gz")
V = FunctionSpace(mesh, "Lagrange", 1)
NU = Function(V)
BB = Function(V)
MU = Function(V)
WW = Function(V)
RHO = Function(V)
File("%s/NU_50.xml"  % data_dir) >> NU
File("%s/BB_50.xml"  % data_dir) >> BB
File("%s/MU_50.xml"  % data_dir) >> MU
File("%s/WW_50.xml"  % data_dir) >> WW
File("%s/RHO_50.xml" % data_dir) >> RHO

# Compute solution
model.parameters.E0 = float(E0)
solution = solver.solve(model, (NU, BB, MU, WW, RHO))

"""
Degenerate Kerr Expressions

"""

import os
from dolfin import *

# Set up directories
#library_dir  = os.path.dirname(os.path.join(os.path.abspath(__file__),'..','geco'))
library_dir  = os.path.dirname(os.path.abspath(__file__))
geometry_dir = os.path.join(library_dir,"..","geco","geometries")
solution_dir = os.path.join("solutions", "degenerate_kerr")

# Set file names for mesh and markers
mesh_prefix = "halfcircle"
R = 50
N = 512
mesh_file = os.path.join(geometry_dir, "%s_%d_%d_mesh.xml.gz" \
                          % (mesh_prefix, R, N))
mark_file = os.path.join(geometry_dir, "%s_%d_%d_markers.xml.gz" \
                          % (mesh_prefix, R, N))
                         
# Create mesh and function spaces
mesh = Mesh(mesh_file)
V = FunctionSpace(mesh, "Lagrange", 1)
markers = MeshFunction("size_t", mesh, mark_file)

# Create spatial coordinates
x = SpatialCoordinate(mesh)
s = x[0]
z = x[1]

# Define Degenerate Kerr Expressions
r = Expression("sqrt(x[0]*x[0] + x[1]*x[1])")
q = Expression("(1.0 + r)*(1.0 + r)", r=r)
nu = Expression("0.5*log( (r*r*q + x[1]*x[1])/( (1.0 + q)*(1.0 + q) - x[0]*x[0] ))", r=r, q=q)
bb = Expression("1.0")
mu = Expression("0.5*log( q/(r*r) + x[1]*x[1]/(r*r*r*r) )", r=r, q=q)
ww = Expression("2.0*(1 + r)/( (1.0 + q)*(1.0 + q) - x[0]*x[0] )", r=r, q=q)

NU = project(nu, V)
BB = project(bb, V)
MU = project(mu, V)
WW = project(ww, V)

# Save metric fields
suffix = "degenerate_kerr"
File(os.path.join(solution_dir, "NU_%d_%s.pvd"  % (R, suffix))) << NU
File(os.path.join(solution_dir, "BB_%d_%s.pvd"  % (R, suffix))) << BB
File(os.path.join(solution_dir, "MU_%d_%s.pvd"  % (R, suffix))) << MU
File(os.path.join(solution_dir, "WW_%d_%s.pvd"  % (R, suffix))) << WW


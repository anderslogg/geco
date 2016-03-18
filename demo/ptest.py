from fenics import *

mesh = UnitSquareMesh(16, 16)
V = FunctionSpace(mesh, "P", 1)
u = Function(V)
A = assemble(TrialFunction(V)*TestFunction(V)*dx)
b = assemble(TestFunction(V)*dx)
x = u.vector()

solve(A, x, b, "gmres")

"""
Equation: Einstein-Vlasov
Ansatz:   EV-E-Polytropic-L-Andreasson (spindle)
          EV-E-Polytropic-L-Polytropic (torus) 

This creates an hoag-type object with a spindle center and a toroidal ring and zero net angular momentum.
The shape is stable under changes to the parameters.
The present choice gives a nice near-vacuum region between the components.
"""

from geco import *

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
solver.parameters.discretization.radius = 100

# Create ansatzes
S = MaterialModel("EV-E-Polytropic-L-Andreasson")
T = MaterialModel("EV-E-Polytropic-L-Polytropic")
model = Constant(0.2)*S + Constant(0.8)*T

# Ansatz parameters
S.parameters.E0 = 0.940
S.parameters.Q  = 3.0
S.parameters.k  = 0.3
S.parameters.l  = 0.5

T.parameters.E0 = 0.940
T.parameters.L0 = 2.4
T.parameters.k = 0.0
T.parameters.l = 0.25

# Set rotation to zero
S.parameters.rotation = False
T.parameters.rotation = False

# Compute solution for initial guess
solution = solver.solve(model)

# Extract solution components
NU, BB, MU, WW, RHO, data = solution

# Plot density
plot(RHO)
interactive()

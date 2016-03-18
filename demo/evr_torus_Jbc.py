"""
Equation: Einstein-Vlasov (rotating)
Ansatz:   EV-E-Polytropic-L-Gaussian

This creates a rotating disk of particles with positive angular momentum.
The angular momentum of the boundary condition is modified at each step of an iteration.
"""

from geco import *
from dolfin import *
from numpy import linspace

# Create solver
solver = EinsteinVlasovSolver()
solver.parameters.output.plot_solution = False
#solver.parameters.discretization.radius = 50
#solver.parameters.discretization.num_steps = 25
#solver.parameters.discretization.resolution = 128

# Create ansatz for initial guess
model = MaterialModel("EV-E-Polytropic-L-Polytropic")
model.parameters.E0 = 0.942
model.parameters.k = 0.0
model.parameters.l = 0.0

# Compute solution for initial guess
solution = solver.solve(model)

# Parameters for J iteration
maxjter = 25
Jtol = 1e-4
totJ = 0.0

# Parameters for increasing L0
L0start = 1.00
L0stop = 2.0
l0vals = linspace(L0start,L0stop, 10)

# Increase L0 value
for l in l0vals:
    print l
    model.parameters.L0 = l

    # Loop over solve modifying J each time.
    for jter in xrange(maxjter):

        info("jter is %d" % jter)

        # Set angular momentum
        info("setting angular momentum to totJ = %g" % totJ)
        solver.parameters.discretization.ang_mom = totJ
        totJp = totJ

        # Solve
        solution = solver.solve(model, solution)

        # Extract solution components
        NU, BB, MU, WW, RHO, data = solution

        # Get angular momentum
        totJ = data['total_angular_momentum']
        info("Computed total J to be %g" % totJ)

        # Check for convergence
        if abs(totJ - totJp)  < Jtol and jter > 0:
            break

# Copy fields to         
fields  = [v.copy(deepcopy=True) for v in solution[0:5]]
                
# Find J = 0 bc solution for comparison
solver.parameters.discretization.ang_mom = 0.0
solution0 = solver.solve(model, solution)
data0 = solution0[6]
totJ0 = data0['total_angular_momentum']


# Check whether iteration converged
if jter == maxjter - 1:
    error("Angular momentum iteration did not converge.")
print
print "Angular momentum iterations converged to within a tolerance of %g." % Jtol
print "Number of iterations was %g." % jter    
    
# Compare fields
fields0 = [v for v in solution0[0:5]]
zero_norms = [sqrt(assemble((u0)**2.0*dx)) for u0 in fields0]
normal_diffs = [sqrt(assemble((u-u0)**2.0*dx))/sqrt(assemble((u0)**2.0*dx)) for u,u0 in zip(fields, fields0)]
field_diffs = [sqrt(assemble((u-u0)**2.0*dx)) for u,u0 in zip(fields, fields0)]

info("Normalized L^2 norm of the difference of the metric fields = %.3g, %.3g, %.3g, %.3g, %.3g" % tuple(normal_diffs))
info("L^2 norm of the difference of the metric fields = %.3g, %.3g, %.3g, %.3g, %.3g" % tuple(field_diffs))
info("L^2 norm of the J=0 BC metric fields = %.3g, %.3g, %.3g, %.3g, %.3g" % tuple(zero_norms))
info("Final angular momentum is %g" % totJ)
info("J = 0 BC angular momentum is %g" % totJ0)    
    
# Plot density
plot(RHO)
interactive()

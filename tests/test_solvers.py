from geco import *

eps = 1e-10


def test_vp():
    "Test VP solver"

    # Create solver
    solver = VlasovPoissonSolver()
    solver.parameters.output.plot_solution = False
    solver.parameters.output.plot_iteration = False
    solver.parameters.discretization.domain_radius = 50
    solver.parameters.discretization.resolution = 16

    # Create ansatz for initial guess
    model = MaterialModel("VP-E-Polytropic-L-Polytropic")

    # Compute solution
    U, RHO, data = solver.solve(model)

    # Check results
    C = data["ansatz_coefficient"]
    assert(abs(C - 0.0025108107357906897) < eps)


def test_ev():
    "Test EV solver, non-rotating"

    # Create solver
    solver = EinsteinVlasovSolver()
    solver.parameters.output.plot_solution = False
    solver.parameters.output.plot_iteration = False
    solver.parameters.discretization.domain_radius = 50
    solver.parameters.discretization.resolution = 16

    # Create ansatz for initial guess
    model = MaterialModel("EV-E-Polytropic-L-Polytropic")
    model.parameters.E0 = 0.942
    model.parameters.k = 0.0
    model.parameters.l = 0.0
    model.parameters.rotation = False

    # Compute solution
    solution = solver.solve(model)

    # Check results
    data = solution.data
    C = data["ansatz_coefficient"]
    assert(abs(C - 0.0007272865499169678) < eps)
    


def test_evr():
    "Test EV solver, rotating"

    # Create solver
    solver = EinsteinVlasovSolver()
    solver.parameters.output.plot_solution = False
    solver.parameters.output.plot_iteration = False
    solver.parameters.discretization.domain_radius = 50
    solver.parameters.discretization.resolution = 16

    # Create ansatz
    model = MaterialModel("EV-E-Polytropic-L-Polytropic")
    model.parameters.E0 = 0.942
    model.parameters.k = 0.0
    model.parameters.l = 0.0

    # Compute solution
    solution = solver.solve(model)

    # Check results
    data = solution.data
    C = data["ansatz_coefficient"]
    assert(abs(C - 0.0015031671809747657) < eps)
    
    
def test_evr_gaussian():
    "Test EV solver, rotating, with E-Poly : L-Guassian ansatz"
    
    # Create solver
    solver = EinsteinVlasovSolver()
    solver.parameters.output.plot_solution = False
    solver.parameters.output.plot_iteration = False    
    solver.parameters.discretization.domain_radius = 50
    solver.parameters.discretization.resolution = 32

    # Create ansatz
    model = MaterialModel("EV-E-Polytropic-L-Gaussian")
    model.parameters.E0 = 0.942
    model.parameters.pm = 1.0
    model.parameters.k  = 0.5
    model.parameters.L0 = 2.0

    # Compute solution
    solution = solver.solve(model)

    # Check results
    data = solution.data
    C = data["ansatz_coefficient"]
    assert(abs(C - 0.015902229428470693) < eps)


if __name__ == "__main__":

    # For debugging individual solvers
    import sys
    tests = {"vp": test_vp, "ev": test_ev, "evr": test_evr, "evr_gaussian": test_evr_gaussian}
    if len(sys.argv) > 1 and sys.argv[1] in tests:
        tests[sys.argv[1]]()

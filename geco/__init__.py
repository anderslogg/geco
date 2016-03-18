from dolfin import *

# Import model classes
from models import MaterialModel
from models import model_data

# Import solver classes
from vpsolver import VlasovPoissonSolver
from evsolver import EinsteinVlasovSolver

# List of solvers (labels)
solvers = ("Vlasov-Poisson", "Einstein-Vlasov")

# Temporary import of temporary class until we make all solvers adaptive
from avpsolver import AdaptiveVlasovPoissonSolver

# List of ansatzes (labels)
ansatzes = model_data

# Mapping from solver label to solver name (prefix)
def solver_name(solver_label):
    return {"Vlasov-Poisson":  "vp",
            "Einstein-Vlasov": "ev"}[solver_label]

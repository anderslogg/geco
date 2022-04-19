from dolfin import *

# Import model classes
from geco.models import MaterialModel
from geco.models import model_data

# Import solver classes
from geco.vpsolver import VlasovPoissonSolver
from geco.evsolver import EinsteinVlasovSolver

# List of solvers (labels)
solvers = ("Vlasov-Poisson", "Einstein-Vlasov")

# Temporary import of temporary class until we make all solvers adaptive
from geco.avpsolver import AdaptiveVlasovPoissonSolver

# List of ansatzes (labels)
ansatzes = model_data

# Mapping from solver label to solver name (prefix)
def solver_name(solver_label):
    return {"Vlasov-Poisson":  "vp",
            "Einstein-Vlasov": "ev"}[solver_label]

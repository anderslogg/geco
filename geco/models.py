"This module defines models for the density."

import os
from dolfin import Expression, error

# List of all available models (grab all files from directory)
library_dir = os.path.dirname(os.path.abspath(__file__))
cppcode_dir = os.path.join(library_dir, "cppcode")
model_data = [m.split(".")[0] for m in os.listdir(cppcode_dir) if \
              m.endswith(".h") and \
             (m.startswith("VP-") or
              m.startswith("EV-"))]

def MaterialModel(model):
    "Create given material model"

    # Fixes subtle issue when model comes in as unicode and should be a string.
    # This pops up when running from the GUI and then calling Parameters.rename.
    model = str(model)

    # Get model data
    if not model in model_data:
        error("Unknown material model: \"%s\"." % str(model))
    model_filename = model + ".h"
    template_filename = model.split("-")[0] + "Ansatz.h"

    # Get library directory
    library_dir = os.path.dirname(os.path.abspath(__file__))
    cppcode_dir = os.path.join(library_dir, "cppcode")

    # Read code
    template_code = open(os.path.join(cppcode_dir, template_filename)).read()
    model_code = open(os.path.join(cppcode_dir, model_filename)).read()

    # Extract relevant model code
    member_functions = model_code.split("// Member functions")[1].split("// Member variables")[0]
    member_variables = model_code.split("// Member variables")[1]

    # Stick specialized code into template and return
    # FIXME: python3 version:
    # template_code.format(member_functions=member_functions, member_functions=member_functions),
    # and in XXAnsatz.h as {member_functions} and {member_variables}
    cppcode = template_code % {"member_functions": member_functions,
                               "member_variables": member_variables}

    # Create expression
    rho = Expression(cppcode=cppcode, degree=1)

    # Set name of ansatz
    rho.parameters.rename(model)

    return rho

# Function for creating 2D representation density (extension to R2 from first quadrant)
def Density2D(rho):

    # Get library directory
    library_dir = os.path.dirname(os.path.abspath(__file__))
    cppcode_dir = os.path.join(library_dir, "cppcode")

    # Read code from file
    cppcode = open(os.path.join(cppcode_dir, "Density2D.h")).read()

    # Build extension module
    rho2d = Expression(cppcode=cppcode, degree=1)

    # Set density
    rho2d.set_density(rho)

    return rho2d

# Function for creating 3D representation density
def Density3D(rho):

    # Get library directory
    library_dir = os.path.dirname(os.path.abspath(__file__))
    cppcode_dir = os.path.join(library_dir, "cppcode")

    # Read code from file
    cppcode = open(os.path.join(cppcode_dir, "Density3D.h")).read()

    # Build extension module
    rho3d = Expression(cppcode=cppcode, degree=1)

    # Set density
    rho3d.set_density(rho)

    return rho3d

# Function for creating indicator function on matter support
def SupportBump(rho):

    # Get library directory
    library_dir = os.path.dirname(os.path.abspath(__file__))
    cppcode_dir = os.path.join(library_dir, "cppcode")

    # Read code from file
    cppcode = open(os.path.join(cppcode_dir, "SupportBump.h")).read()

    # Build extension module
    rho_support = Expression(cppcode=cppcode, degree=1)

    # Set density
    rho_support.set_density(rho)

    return rho_support

# Function for creating point cloud representation of density
def PointCloud(rho, R, M, resolution, num_points):

    # Get library directory
    library_dir = os.path.dirname(os.path.abspath(__file__))
    cppcode_dir = os.path.join(library_dir, "cppcode")

    # Read code from file
    cppcode = open(os.path.join(cppcode_dir, "PointCloud.h")).read()

    # Build extension module
    point_cloud = Expression(cppcode=cppcode, degree=1)

    # Set density
    point_cloud.set_parameters(rho, R, M, resolution, num_points)

    return point_cloud

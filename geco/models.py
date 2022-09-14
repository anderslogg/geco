# Copyright 2019 Anders Logg, Ellery Ames, Haakan Andreasson
#
# This file is part of GECo. GECo is free software: you can
# redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# GECo is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GECo. If not, see <https://www.gnu.org/licenses/>.

"""
This module defines models for the density.
"""

import os

from dolfin import CompiledExpression, compile_cpp_code, error

# List of all available models (grab all files from directory)
library_dir = os.path.dirname(os.path.abspath(__file__))
cppcode_dir = os.path.join(library_dir, "cppcode")
model_data = [
    m.split(".")[0]
    for m in os.listdir(cppcode_dir)
    if m.endswith(".h") and (m.startswith("VP-") or m.startswith("EV-"))
]

# Cache for generated Ansatz classes
model_cache = {}


def MaterialModel(model):
    "Create given material model"

    # Fixes subtle issue when model comes in as unicode and should be a string.
    # This pops up when running from the GUI and then calling Parameters.rename.
    model = str(model)

    # Get model data
    if model not in model_data:
        error('Unknown material model: "%s".' % str(model))
    
    # Select appropriate bindings file and module name
    bindings_filename = model.split('-')[0] + "Bindings.h"
    model_short_name = "".join(model.split('-'))
    module_name = "compiled_module." + model_short_name

    # Get library directory
    library_dir = os.path.dirname(os.path.abspath(__file__))
    cppcode_dir = os.path.join(library_dir, "cppcode")

    # Read code
    cppcode = open(os.path.join(cppcode_dir, bindings_filename)).read()

    # Get compiled Anzatz class
    if model in model_cache:
        print('Reusing class %s from cache' % model)
        compiled_model = model_cache[model]
    else:
        print('Compiling new class %s' % model)
        compiled_module = compile_cpp_code(cppcode, include_dirs=[cppcode_dir])
        compiled_model = eval(module_name)
        model_cache[model] = compiled_model

    # Create compiled expression
    rho = CompiledExpression(compiled_model(), degree=1)

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

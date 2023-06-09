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
A library of functions to compute various characteristics of the
numerically computed solutions.
"""

import os

import numpy as np
from dolfin import *

# Default quantities called by evsolver
default_quantities = [
    "total_mass",
    "rest_mass",
    "total_angular_momentum",
    "fractional_binding_energy",
    "ergo_region",
]


def compute_default_physical_quantities(U):
    return compute_physical_quantities(default_quantities, U)


def compute_physical_quantities(quantities, U):
    return dict((q, compute_physical_quantity(q, U)) for q in quantities)


def compute_physical_quantity(quantity, U):
    return globals()["compute_" + quantity](U)


# Total mass
def compute_total_mass(U):

    x = SpatialCoordinate(U.mesh)
    r = x[0]
    mass = assemble(2 * 2 * pi * U.RHO * r * dx(U.mesh))

    return mass


# Rest mass
def compute_rest_mass(U):

    x = SpatialCoordinate(U.mesh)
    r = x[0]
    rest_mass = assemble(2 * 2 * pi * U.RMD * r * dx(U.mesh))

    return rest_mass


# Total angular momentum
def compute_total_angular_momentum(U):

    x = SpatialCoordinate(U.mesh)
    r = x[0]
    J = assemble(
        -2
        * 2
        * pi
        * exp(-4.0 * U.NU)
        * U.BB
        * (U.P03 + (r * U.BB) ** 2 * U.WW * U.P33)
        * r
        * dx(U.mesh)
    )

    return J


# The tt-component of the metric tensor
def compute_gtt_metric_component(U):

    x = SpatialCoordinate(U.mesh)
    r = x[0]
    gtt = project(
        -exp(2.0 * U.NU) * (1.0 - (U.WW * r * U.BB) ** 2 * exp(-4.0 * U.NU)), U.V
    )

    return gtt


# The maximum of gtt
def compute_gtt_max(U):
    gtt = compute_gtt_metric_component(U)
    gtt_max = gtt.vector().max()

    return gtt_max


# Test for ergo region
def compute_ergo_region(U):

    gtt = compute_gtt_metric_component(U)
    gtt_max = gtt.vector().max()

    return gtt_max > 0


# Compute fractional binding energy
def compute_fractional_binding_energy(U):

    #    m  = compute_total_mass(U)
    #    rm = compute_rest_mass(U)

    try:
        m = float(U.data["prescribed_mass"])
    except KeyError:
        m = float(U.data["mass"])
    except:
        print("No key for mass found")
        raise

    try:
        rm = float(U.data["rest_mass"])
    except:
        rm = compute_rest_mass(U)

    return 1.0 - m / rm


# Reflection plane radii of support
def compute_reflection_plane_support(U):

    r_supp = float(U.data["radius_of_support"])

    # compute r-values
    r_res = 10000
    thresh = 1e-5
    rvals = np.linspace(0, r_supp, r_res)
    RHOvals = np.array([U.RHO(r, 0) for r in rvals])
    support = np.where(RHOvals > thresh)[0]
    try:
        r_inner = rvals[min(support)]
        r_outer = rvals[max(support)]
        r_peak = rvals[np.argmax(RHOvals)]
    except ValueError:
        # print ("Issue with matter support, setting NA value for rho.")
        if min(RHOvals) < thresh:
            print("Density is below threshold, setting rvals to 0.0")
        r_inner = 0.0
        r_outer = 0.0
        r_peak = 0.0

    return [r_inner, r_peak, r_outer]


def compute_areal_radius_of_support(U, r_supp):

    m = compute_total_mass(U)
    R0 = r_supp * (1.0 + m / (2.0 * r_supp)) ** 2

    return R0


# Rcirc, the length of the axisymmetric Killing vector field
def compute_Rcirc_func(U):
    x = SpatialCoordinate(U.mesh)
    r = x[0]
    RcircF = project(r * U.BB * exp(-U.NU), U.V)

    return RcircF


def compute_Rcirc_values(U, rvals):

    # Compute the function
    RcircF = compute_Rcirc_func(U)

    return [RcircF(r, 0.0) for r in rvals]


# Mass aspect function
def compute_mass_aspect(U):

    r_supp = float(U.data["radius_of_support"])

    Rcirc_func = compute_Rcirc_func(U)

    x = SpatialCoordinate(U.mesh)
    r = x[0]

    sres = 100
    slist = np.linspace(0.0, r_supp, sres)
    masslist = []

    for s in slist:
        g = Expression("s - x[0]", degree=1, s=s)
        eta = (
            2 * 2 * pi * assemble(conditional(ge(g, 0), 1, 0) * U.RHO * r * dx(U.mesh))
        )  # Doubled to account for quarter disk
        eta *= 2  # 2m/r, is less than 1 for regular solutions
        eta /= Rcirc_func(s, 0.0)
        masslist.append([s, eta])

    maximum = max(masslist, key=lambda x: x[1])

    return maximum


# Lapse
def compute_lapse_func(U):
    lapse = project(exp(U.NU), U.V)

    return lapse


def compute_lapse_values(U, rvals):

    # Get domain radius
    try:
        r_boundary = U.data["domain_radius"]
    except KeyError:
        r_boundary = U.data["radius"]
    try:
        r_boundary = r_boundary[-1]
    except IndexError:
        pass

    lapse = compute_lapse_func(U)

    # Set any non-float rvals to the boundary
    for r in rvals:
        if not isinstance(r, float):
            r = float(r_boundary)

    return [lapse(r, 0.0) for r in rvals]


# Scaled ZAMO redshift
## Function
def compute_sZAMO_redshift_func(U):
    szRS_func = project(1.0 - exp(U.NU), U.V)

    return szRS_func


## Evaluated at points in the reflection plane
def compute_sZAMO_redshift_values(U, rvals):

    # Get domain radius
    try:
        r_boundary = U.data["domain_radius"]
    except KeyError:
        r_boundary = U.data["radius"]
    try:
        r_boundary = r_boundary[-1]
    except IndexError:
        pass

    zamo_RS_func = compute_sZAMO_redshift_func(U)

    # Set any non-float rvals to the boundary
    for r in rvals:
        if not isinstance(r, float):
            r = float(r_boundary)

    return [zamo_RS_func(r, 0.0) for r in rvals]


# Linear energy density
def compute_linear_energy_density(U):

    led = project(
        U.P00 + U.WW * exp(-4 * U.NU) * U.P03, U.V
    )  # -T_t^t mult. by exp(2MU) for the integration measure.
    linearED = 2 * assemble(led * dx(U.mesh))  # mult. times 2 b/c of quarter disk

    return linearED


# Azimuthal pressure
def compute_azimuthal_pressure(U):

    x = SpatialCoordinate(U.mesh)
    r = x[0]

    tmp = project(
        -U.WW * exp(-4 * U.NU) * U.P03
        + (1.0 - (r * U.BB * U.WW) ** 2 * exp(-4 * U.NU)) * U.P33,
        U.V,
    )  # T_phi^phi mult. by exp(2MU) for the integration measure.
    azimuthal_pressure = 2 * assemble(
        tmp * dx(U.mesh)
    )  # mult. times 2 b/c of quarter disk

    return azimuthal_pressure


# Meridional plane pressure
def compute_meridional_pressure(U):

    # x = SpatialCoordinate(U.mesh)
    # r = x[0]

    # tmp = project( U.P11 , U.V) # T_rho^rho + T_z^z mult. by exp(2MU) for the integration measure.
    meridional_pressure = 2 * assemble(
        U.P11 * dx(U.mesh)
    )  # mult. times 2 b/c of quarter disk

    return meridional_pressure


# Hakan's parameter
def compute_hakan_parameter(U):

    r_peak = compute_reflection_plane_support(U)[1]

    return exp(2 * U.NU(r_peak, 0.0)) / U.BB(r_peak, 0.0)


# Dain Mass
def compute_dain_mass(U):

    x = SpatialCoordinate(U.mesh)
    r = x[0]
    matter_terms = (
        2 * pi * (U.P00 - U.WW * U.WW * (r * U.BB) ** 2 * exp(-4.0 * U.NU) * U.P33)
    )
    field_terms = (
        dot(grad(U.BB), grad(U.BB)) / (4.0 * U.BB * U.BB)
        + dot(grad(U.NU), grad(U.NU)) / 4.0
        - dot(grad(U.BB), grad(U.NU)) / (2.0 * U.BB)
        + (r * U.BB / 4) ** 2 * exp(-4.0 * U.NU) * dot(grad(U.WW), grad(U.WW))
    )
    dain_mass = assemble((matter_terms + field_terms) * r * dx(U.mesh))

    return dain_mass

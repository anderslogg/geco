#!/usr/bin/env python
#
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
A collection of tools for plotting data from GECo and adaptivesolver.py
solutions or a sequence of such solutions.

Usage:
------
list_data(data_files)
 # print list of data that can be plotted

list_derived_data()
 # print list of derived quantities available to plot

gecoplot(data_files, 'E0', 'frac_binding_energy', labels='ergo_region', legend=None, converged_only=True, verbose=False)
 # plots Eb vs E0 and labels the points by whether they contain an ergo_region.
 # only plots converged solutions.
 # won't warn when data is missing from certain files. To change set verbose=True
"""

import glob
import os
import sys

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams

axes_label_dict = {
    "E0": "$E_0$",
    "L0": "$L_0$",
    "Rcirc_inner": "$\bar R_{circ}$ at inner radius",
    "Rcirc_peak": "$R_{circ}$ at radius of peak",
    "Rcirc_outer": "$\overline{R}_{circ}$ at outer radius",
    "Rcirc_max": r"$R_{circ}$ at maximum of $\Gamma$",
    "adaptive_theta": "Adaptive Damping Parameter",
    "anderson_depth": "Anderson Depth",
    "ansatz_coefficient": "Ansatz Coefficient",
    "areal_radius_of_support": "Radius of support in areal-like coordinates",
    "central_redshift": "Central Redshift $Z_c$",
    "degree": "Degree of the ???",
    "ergo_region": "Indicator of an ergo region",
    "fractional_binding_energy": "Fractional Binding Energy $E_b$",
    "gamma": "A compactness parameter $2M/R$",
    "gtt_max": "Maximum of the $g_{tt}$ component",
    "k": "Exponent of the energy part of the distribution",
    "l": "Exponent of the momentum part of the distribution",
    "krylov_tolerance": "Krylov Tolerance",
    "mass": "Total Mass of the Solution",
    "zamo_redshift_outer": r"$\bar{Z}_o$",
    "zamo_redshift_inner": r"$\bar{Z}_i$",
    "zamo_redshift_peak": r"$\bar{Z}_p$",
    "zamo_redshift_max": r"$\bar{Z}_\Gamma$",
    "zamo_redshift_orig": r"$\bar{Z}_{origin}$",
    "particle_mass": "Particle Mass",
    "r_inner": "Inner radius of distribution",
    "r_outer": "Outer radius of distribution",
    "r_peak": "Radius of peak of distribution",
    "radius_of_support": "Coordinate Radius of Support",
    "rest_mass": "Rest Mass",
    "solution_converged": "Indicate the solution converged",
    "total_angular_momentum": "Total Angular Momentum",
    "ri/ro": "Inner radius of support over outer radius of support",
    "M_squared_over_J": "$M^2 / J$",
    "J_over_M_squared": "$J/M^2$",
    "M_over_Rcirc": "Mass over Rcirc",
    "mass_aspect_max": r"$\Gamma =$ max($2m/R_{circ}(z=0)$)",
    "mass_aspect_max_r": "Radius of maximum of $2m/R_{circ}$",
    "central_lapse": "Central Lapse",
    "peak_lapse": "Lapse at matter peak",
    "Rcirc_squared_over_J": "Rcirc squared over J",
    "HMV_deficit_angle": "$4\pi (u + m - k)$",
    "meridional_pressure": "Meridional Pressure",
    "fourpi_meridional_pressure": "$4\pi m$",
}


derived_quantities = {
    "ri/ro": [["r_inner", "r_outer"], "df_radius_ratio"],
    "normalized_central_redshift": [
        ["central_redshift"],
        "df_normalized_central_redshift",
    ],
    "M_squared_over_J": [["mass", "total_angular_momentum"], "df_M_squared_over_J"],
    "J_over_M_squared": [["mass", "total_angular_momentum"], "df_J_over_M_squared"],
    "M_over_Rcirc": [["mass", "Rcirc"], "df_M_over_Rcirc"],
    "Rcirc_squared_over_J": [
        ["Rcirc", "total_angular_momentum"],
        "df_Rcirc_squared_over_J",
    ],
    "HMV_deficit_angle": [
        ["linear_energy_density", "azimuthal_pressure", "meridional_pressure"],
        "df_deficit_angle_estimate",
    ],
    "fourpi_meridional_pressure": [
        ["meridional_pressure"],
        "df_fourpi_meridional_pressure",
    ],
}


save_dir = os.path.dirname(os.path.realpath(__file__))


def get_data_index(data_file, data):
    "Returns the column number corresponding to data in the solution data_file"

    header = np.genfromtxt(data_file, max_rows=1, dtype=str, delimiter=",").tolist()
    #    print(header)
    #    names = header[0,:].tolist()
    #    names = [s[2:-1] for s in names] # strips b'...' in python3
    try:
        data_index = header.index(data)
        return data_index
    except ValueError:
        raise ValueError


def list_data(data_files):
    # Lists data available for plotting

    flat_list = [item for sublist in data_files for item in sublist]

    headers = []
    warning = False

    # Find union of all data quantities
    for data_file in flat_list:
        header = np.genfromtxt(data_file, max_rows=1, delimiter=",", dtype=str)
        headers = list(set(headers).union(header))

        # If non-overlapping, print warning.
        if not len(headers) == len(header):
            warning = True
        elif not all(header == sorted(headers)):
            warning = True

    for h in np.sort(headers):
        print(h)

    if warning:
        print()
        print("Warning: Not all data files contain all the listed data")


def list_derived_data():
    # Lists derived data available for plotting

    # Find union of all data quantities
    for k in derived_quantities.keys():
        print(k)


def look_up_labels(xdata, ydata, labels):
    # Returns names for xdata, ydata, and labels

    label_name = ""
    xlabel = ""
    ylabel = ""
    try:
        xlabel = axes_label_dict[xdata]
        ylabel = axes_label_dict[ydata]

        if labels is not None:
            label_name = axes_label_dict[labels]
    except KeyError:
        print("Names for requested quantity have not been entered in  axes_label_dict.")
        print("Returning blank names for some requested quantities.")

    return xlabel, ylabel, label_name


def set_save_dir(save_directory):
    # Sets directory in which to store files
    global save_dir
    save_dir = save_directory
    # print('files will be saved to %s' %save_dir)


def get_data(data_file, data_name):
    # returns data from file if available

    # retrieve names of available data
    header = np.genfromtxt(data_file, max_rows=1, delimiter=",", dtype=str)

    # if requested data is in list of names retrieve the data
    if data_name in header:

        if data_name == "ergo_region" or data_name == "solution_converged":
            typ = bool
        else:
            typ = float

        data_index = get_data_index(data_file, data_name)
        data = np.genfromtxt(
            data_file,
            delimiter=",",
            skip_header=1,
            dtype=typ,
            usecols=(data_index),
            unpack=True,
        )

    # else look in the list of available derived quantities
    elif data_name in derived_quantities.keys():

        # print('  Retrieving data for derived quantity...')

        # get fundamental data
        fdata = derived_quantities[data_name][0]
        derived_data = []
        for dname in fdata:
            data_index = get_data_index(data_file, dname)
            data = np.genfromtxt(
                data_file,
                delimiter=",",
                skip_header=1,
                usecols=(data_index),
                unpack=True,
            )
            data = np.array([data])
            derived_data = np.append(derived_data, data)

        derived_data = derived_data.reshape(len(fdata), len(data))

        # compute derived thing according to some formula... (also stored in derived _dict?)
        func_handle = derived_quantities[data_name][1]
        if func_handle == "df_radius_ratio":
            data = df_radius_ratio(derived_data)
        elif func_handle == "df_normalized_central_redshift":
            data = df_normalized_central_redshift(derived_data)
        elif func_handle == "df_M_squared_over_J":
            data = df_M_squared_over_J(derived_data)
        elif func_handle == "df_J_over_M_squared":
            data = df_J_over_M_squared(derived_data)
        elif func_handle == "df_M_over_Rcirc":
            data = df_M_over_Rcirc(derived_data)
        elif func_handle == "df_Rcirc_squared_over_J":
            data = df_Rcirc_squared_over_J(derived_data)
        elif func_handle == "df_deficit_angle_estimate":
            data = df_deficit_angle_estimate(derived_data)
        elif func_handle == "df_fourpi_meridional_pressure":
            data = df_fourpi_meridional_pressure(derived_data)
        else:
            print("A function for this quantity has not been defined.")

    # else print an warning message
    else:
        print("Quantity not found. Returning zero")
        data = [0.0]

    return np.array(data).reshape(-1)


def geco_pp_plot(
    data_runs,
    xdata_name,
    ydata_name,
    legend_labels=None,
    point_labels=None,
    markers=None,
    ylims=None,
    savefig=False,
):
    # Specialized for plotting postprocessing data in that each data file should contain only one line of data
    # (no unconverged solutions present)
    # Takes a list of data runs data_runs = [ [run1_file1, run1_file2, ...], [run2_file1, run2_file2, ...], ... ]

    # Set figure size
    my_dpi = 300
    plt.figure(figsize=(8, 4), dpi=my_dpi)
    #    matplotlib.rcParams['figure.figsize'] = (2.0, 1.0)
    #    legend_params = {'legend.fontsize':4, 'axes.labelsize':6}
    #    plt.rcParams.update(legend_params)
    #    plt.rc('text', usetex=True)

    # Set legend labels
    if legend_labels == None:
        run_labels = ["Run {:}".format(i) for i in range(len(data_runs))]
    else:
        run_labels = legend_labels

    # Set markers
    if markers == None:
        run_markers = ["o" for i in range(len(data_runs))]
    else:
        run_markers = markers

    # Generate x and y data in each run
    for data_run, run_label, run_marker in zip(data_runs, run_labels, run_markers):

        # Generate xdata and ydata
        xdata = np.array([get_data(f, xdata_name)[0] for f in data_run])
        ydata = np.array([get_data(f, ydata_name)[0] for f in data_run])

        # Make sure data corresponds to converged solution
        conv_data = np.array([get_data(f, "solution_converged")[0] for f in data_run])
        xdata = xdata[np.where(conv_data == True)]
        ydata = ydata[np.where(conv_data == True)]

        # plot data
        plt.plot(
            xdata,
            ydata,
            marker=run_marker,
            markersize=4,
            linestyle="None",
            label=run_label,
        )

    # Look up axes point_labels
    xlabel, ylabel, label_name = look_up_labels(xdata_name, ydata_name, point_labels)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.tick_params(tickdir="in", length=2, width=2, labelsize=16)
    plt.grid()
    if not legend_labels == "empty":
        plt.legend(fontsize=16)

    if ylims != None:
        try:
            ybot, ytop = ylims
            plt.ylim(ybot, ytop)
        except:
            print("Give (ybot, ytop)")
            raise ValueError

    # Save file if desired
    if savefig:
        file_name = "%s_vs_%s.png" % (ydata_name, xdata_name)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        save_file = os.path.join(save_dir, file_name)
        print("Saving figure as %s" % save_file)
        plt.savefig(save_file, dpi=my_dpi, bbox_inches="tight")

    plt.show()
    # return plt


def gecoplot(
    data_runs,
    xdata,
    ydata,
    point_labels=None,
    converged_only=True,
    savefig=False,
    verbose=False,
):
    # plots ydata vs xdata for data in data_files.
    # Capable of plotting any unconverged data.
    # Options: labels, legend, converged_only, savefig
    # TODO: add legend and different coloring abilities. Might require more structure in the input files though...

    # Set figure size
    my_dpi = 96
    plt.figure(figsize=(1152 / my_dpi, 576 / my_dpi), dpi=my_dpi)
    # matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)

    for data_run in data_runs:

        run_xdata = list()
        run_ydata = list()

        for data_file in data_run:

            # look up index for requested data
            try:

                x_data = get_data(data_file, xdata)
                y_data = get_data(data_file, ydata)

                if (
                    point_labels == "ergo_region"
                    or point_labels == "solution_converged"
                ):
                    label_data = get_data(data_file, point_labels)

                elif point_labels is not None:
                    label_data = np.round(get_data(data_file, point_labels), 3)

                elif point_labels is None:
                    label_data = ""

                if converged_only:
                    sc_data = get_data(data_file, "solution_converged")

            except ValueError:
                if verbose:
                    print("Desired data not found in '%s'" % (data_file))
                continue

            if converged_only:
                # drop un-converged solutions
                x_data = x_data[np.where(sc_data == True)]
                y_data = y_data[np.where(sc_data == True)]

            # FIXME: Point labels are likely broken.
            if point_labels is not None:
                for label, x, y in zip(label_data, x_data, y_data):
                    plt.annotate(
                        str(label),
                        xy=(x, y),
                        xytext=(-2, 2),
                        textcoords="offset points",
                        ha="left",
                        va="bottom",
                    )

            # pack points in to one list
            run_xdata.append(x_data)
            run_ydata.append(y_data)

        # flatten and plot run data
        run_xdata = [x for xa in run_xdata for x in xa]
        run_ydata = [y for ya in run_ydata for y in ya]
        plt.plot(run_xdata, run_ydata, "o", markersize=4)

    # Look up axes point_labels
    xlabel, ylabel, label_name = look_up_labels(xdata, ydata, point_labels)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.grid()

    # Save file if desired
    if savefig:
        file_name = "%s_vs_%s.png" % (ydata, xdata)
        save_file = os.path.join(save_dir, file_name)
        print("Saving figure as %s" % save_file)
        plt.savefig(save_file, dpi=100, bbox_inches="tight")

    plt.show()


# Highlight Data Point
def highlight_point(ax, data_file, xdata_name, ydata_name, hmarker):

    # Get data
    xdata = get_data(data_file, xdata_name)[0]
    ydata = get_data(data_file, ydata_name)[0]

    ax.plot(xdata, ydata, marker=hmarker, markersize=10)

    return ax


# Interesting ratios of data
def df_radius_ratio(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [r_inner, r_outer]
    # where r_inner and r_outer are possibly lists of values.

    r_inner = np.array(arg_array[0]).reshape(-1)
    r_outer = np.array(arg_array[1]).reshape(-1)

    return [ri / ro for ri, ro in zip(r_inner, r_outer)]


def df_normalized_central_redshift(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [central_redshift]

    zc_data = np.array(arg_array[0]).reshape(-1)

    return [zc / (1.0 + zc) for zc in zc_data]


def df_M_squared_over_J(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [Mass, total_angular_momentum]

    mass = np.array(arg_array[0]).reshape(-1)
    angm = np.array(arg_array[1]).reshape(-1)

    return [m * m / j for m, j in zip(mass, angm)]


def df_J_over_M_squared(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [Mass, total_angular_momentum]

    mass = np.array(arg_array[0]).reshape(-1)
    angm = np.array(arg_array[1]).reshape(-1)

    return [j / (m * m) for m, j in zip(mass, angm)]


def df_M_over_Rcirc(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [Mass, Rcirc]

    mass = np.array(arg_array[0]).reshape(-1)
    rcirc = np.array(arg_array[1]).reshape(-1)

    return [m / r for m, r in zip(mass, rcirc)]


def df_Rcirc_squared_over_J(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [Rcirc, J]

    rcirc = np.array(arg_array[0]).reshape(-1)
    angm = np.array(arg_array[1]).reshape(-1)

    return [r * r / j for r, j in zip(rcirc, angm)]


def df_deficit_angle_estimate(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [linear_energy_density, azimuthal_pressure, meridional_pressure]

    linearED = np.array(arg_array[0]).reshape(-1)
    aziPress = np.array(arg_array[1]).reshape(-1)
    merPress = np.array(arg_array[2]).reshape(-1)

    return [4 * np.pi * (l + m - a) for l, a, m in zip(linearED, aziPress, merPress)]


def df_fourpi_meridional_pressure(arg_array):
    # Takes an arg_array consisting of
    # arg_array = [meridional_pressure]

    merPress = np.array(arg_array[0]).reshape(-1)

    return [4 * np.pi * m for m in merPress]


########################################################################
########################################################################
# OLD FUNCTIONS


def normalized_redshift_vs_radius_ratio(sol_dirs, labels):
    "Plots the normalized redshift vs the ri/ro, as in Ansorg et al."

    for sol_dir in sol_dirs:

        try:
            label_index = get_data_index(labels, sol_dir)
            rinner_index = get_data_index("r_inner", sol_dir)
            router_index = get_data_index("r_outer", sol_dir)
            redshift_index = get_data_index("central_redshift", sol_dir)

            label_data, ri_data, ro_data, rs_data = np.loadtxt(
                sol_dir + "/data.csv",
                delimiter=",",
                skiprows=1,
                usecols=(label_index, rinner_index, router_index, redshift_index),
                unpack=True,
            )

        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        la_iter = np.array([label_data]).reshape(-1)
        ri_iter = np.array([ri_data]).reshape(-1)
        ro_iter = np.array([ro_data]).reshape(-1)
        rs_iter = np.array([rs_data]).reshape(-1)
        r_ratio = [rid / rod for rid, rod in zip(ri_iter, ro_iter)]
        zc_norm = [rsd / (1.0 + rsd) for rsd in rs_iter]
        plt.plot(r_ratio, zc_norm, ":ro")
        if labels:
            for label, x, y in zip(la_iter, r_ratio, zc_norm):
                plt.annotate(
                    str(label),
                    xy=(x, y),
                    xytext=(-2, 2),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                )

    plt.xlabel("Inner radius/Outer radius")
    plt.ylabel("Normalized central redshift, $Z_C/(1 + Z_C)$")

    plt.show()


def angular_momentum_vs_E0(sol_dirs):
    "Plots the angular momentum vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            J_index = get_data_index("total_angular_momentum", sol_dir)
        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        e0_data, j_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(e0_index, J_index),
            unpack=True,
        )

        plt.plot(e0_data, j_data, ":ro")

    plt.xlabel("$E_0$")
    plt.ylabel("Total Angular Momentum, $J$")

    plt.show()


def frac_binding_energy_vs_E0(sol_dirs):
    "Plots the fractional binding energy vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            Eb_index = get_data_index("frac_binding_energy", sol_dir)
        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        e0_data, eb_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(e0_index, Eb_index),
            unpack=True,
        )

        plt.plot(e0_data, eb_data, ":ro")

    plt.xlabel("$E_0$")
    plt.ylabel("Fractional Binding Energy, $E_b$")

    plt.show()


def normalized_central_redshift_vs_E0(sol_dirs):
    "Plots the normalized central redshift vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            Zc_index = get_data_index("central_redshift", sol_dir)
        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        e0_data, zc_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(e0_index, Zc_index),
            unpack=True,
        )

        zc_norm = zc_data / (1.0 + zc_data)

        plt.plot(e0_data, zc_norm, ":ro")

    plt.xlabel("$E_0$")
    plt.ylabel("Normalized Central Redshift, $Z_c/(1 + Z_c)$")

    plt.show()


def normalized_central_redshift_vs_frac_binding_energy(sol_dirs, labels):
    "Plots the normalized redshift vs fractional binding energy"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            Zc_index = get_data_index("central_redshift", sol_dir)
            Eb_index = get_data_index("frac_binding_energy", sol_dir)
        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        e0_data, zc_data, eb_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(e0_index, Zc_index, Eb_index),
            unpack=True,
        )

        e0_iter = np.array([e0_data]).reshape(-1)
        eb_iter = np.array([eb_data]).reshape(-1)
        zc_iter = np.array([zc_data]).reshape(-1)
        zc_norm = [rsd / (1.0 + rsd) for rsd in zc_iter]

        plt.plot(eb_iter, zc_norm, ":ro")

        if labels:
            for label, x, y in zip(e0_iter, eb_iter, zc_norm):
                plt.annotate(
                    str(label),
                    xy=(x, y),
                    xytext=(-2, 2),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                )

    plt.xlabel("Fractional Binding Energy $E_b$")
    plt.ylabel("Normalized Central Redshift, $Z_c/(1 + Z_c)$")

    plt.show()


def frac_binding_energy_vs_normalized_central_redshift(sol_dirs, labels):
    "Plots the normalized redshift vs fractional binding energy"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            Zc_index = get_data_index("central_redshift", sol_dir)
            Eb_index = get_data_index("frac_binding_energy", sol_dir)
        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        e0_data, zc_data, eb_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(e0_index, Zc_index, Eb_index),
            unpack=True,
        )

        e0_iter = np.array([e0_data]).reshape(-1)
        eb_iter = np.array([eb_data]).reshape(-1)
        zc_iter = np.array([zc_data]).reshape(-1)
        zc_norm = [rsd / (1.0 + rsd) for rsd in zc_iter]

        plt.plot(zc_norm, eb_iter, ":ro")

        if labels:
            for label, x, y in zip(e0_iter, zc_norm, eb_iter):
                plt.annotate(
                    str(label),
                    xy=(x, y),
                    xytext=(-2, 2),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                )

    plt.ylabel("Fractional Binding Energy $E_b$")
    plt.xlabel("Normalized Central Redshift, $Z_c/(1 + Z_c)$")

    plt.show()


def frac_binding_energy_vs_central_redshift(sol_dirs, labels, converged_only):
    "Plots the normalized redshift vs fractional binding energy"

    for sol_dir in sol_dirs:

        try:
            Ergo_index = get_data_index("ergo_region", sol_dir)
            Zc_index = get_data_index("central_redshift", sol_dir)
            Eb_index = get_data_index("frac_binding_energy", sol_dir)
            Sc_index = get_data_index("solution_converged", sol_dir)

            zc_data, eb_data = np.loadtxt(
                sol_dir + "/data.csv",
                delimiter=",",
                skiprows=1,
                usecols=(Zc_index, Eb_index),
                unpack=True,
                converters={Ergo_index: lambda s: int(s == "True")},
            )

            sc_data = np.genfromtxt(
                sol_dir + "/data.csv",
                delimiter=",",
                skip_header=1,
                usecols=(Sc_index),
                unpack=True,
                dtype=bool,
            )

            if labels is not None:
                La_index = get_data_index(labels, sol_dir)
                if labels == "ergo_region":
                    typ = bool
                else:
                    typ = float
                la_data = np.genfromtxt(
                    sol_dir + "/data.csv",
                    delimiter=",",
                    skip_header=1,
                    usecols=(La_index,),
                    unpack=True,
                    dtype=typ,
                )

        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        eb_iter = np.array([eb_data]).reshape(-1)
        zc_iter = np.array([zc_data]).reshape(-1)
        sc_iter = np.array([sc_data]).reshape(-1)

        if converged_only:
            # drop un-converged solutions
            # converged_data = sc_iter == True
            eb_iter = eb_iter[np.where(sc_iter == True)]
            zc_iter = zc_iter[np.where(sc_iter == True)]

        plt.plot(zc_iter, eb_iter, ":ro")

        if labels is not None:
            la_iter = np.array([la_data]).reshape(-1)
            for label, x, y in zip(la_iter, zc_iter, eb_iter):
                plt.annotate(
                    str(label),
                    xy=(x, y),
                    xytext=(-2, 2),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                )

    plt.ylabel("Fractional Binding Energy $E_b$")
    plt.xlabel("Central Redshift, $Z_c$")
    plt.figure(figsize=(50, 50))

    plt.show()


def normalized_central_redshift_vs_M_squared_over_J(sol_dirs):
    "Plots the normalized redshift vs M^2/J"

    for sol_dir in sol_dirs:

        try:
            Zc_index = get_data_index("central_redshift", sol_dir)
            Jt_index = get_data_index("total_angular_momentum", sol_dir)
        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        zc_data, jt_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(Zc_index, Jt_index),
            unpack=True,
        )

        zc_norm = zc_data / (1.0 + zc_data)
        m_over_j = 1 / (jt_data)

        plt.plot(m_over_j, zc_norm, ":ro")

    plt.xlabel("Mass square over Total Angular Momentum  $M^2/J$")
    plt.ylabel("Normalized Central Redshift, $Z_c/(1 + Z_c)$")

    plt.show()


def M_over_Rcirc_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            Rcirc_index = get_data_index("Rcirc", sol_dir)
            m_index = get_data_index("mass", sol_dir)
        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        e0_data, rc_data, m_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(e0_index, Rcirc_index, m_index),
            unpack=True,
        )

        m_over_rc = m_data / rc_data
        plt.plot(e0_data, m_over_rc, ":ro")

    plt.xlabel("$E_0$")
    plt.ylabel("$M/R_{circ}$")

    plt.show()


def r_peak_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            rpeak_index = get_data_index("r_peak", sol_dir)

            e0_data, rp_data = np.loadtxt(
                sol_dir + "/data.csv",
                delimiter=",",
                skiprows=1,
                usecols=(e0_index, rpeak_index),
                unpack=True,
            )

        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        plt.plot(e0_data, rp_data, ":ro")

    plt.xlabel("$E_0$")
    plt.ylabel("Radius peak density")

    plt.show()


def r_values_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"

    for sol_dir in sol_dirs:

        try:
            e0_index = get_data_index("E0", sol_dir)
            rp_index = get_data_index("r_peak", sol_dir)
            ri_index = get_data_index("r_inner", sol_dir)
            ro_index = get_data_index("r_outer", sol_dir)

            e0_data, rp_data, ri_data, ro_data = np.loadtxt(
                sol_dir + "/data.csv",
                delimiter=",",
                skiprows=1,
                usecols=(e0_index, rp_index, ri_index, ro_index),
                unpack=True,
            )

        except ValueError:
            print("Desired data not found in '%s'" % (sol_dir))
            continue

        plt.plot(e0_data, rp_data, ":ro")
        plt.plot(e0_data, ri_data, ":bo")
        plt.plot(e0_data, ro_data, ":go")

    plt.xlabel("$E_0$")
    plt.ylabel("Peak, Outer, and Inner radius of matter support")
    plt.legend(["peak rad.", "inner rad.", "outer rad."], loc="upper left")

    plt.show()


def ergo_vs_E0(sol_dirs):
    "Plots M/Rcirc vs E0"

    for sol_dir in sol_dirs:

        #  try:
        e0_index = get_data_index("E0", sol_dir)
        er_index = get_data_index("ergo_region", sol_dir)

        e0_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(e0_index,),
            unpack=True,
        )

        er_data = np.loadtxt(
            sol_dir + "/data.csv",
            delimiter=",",
            skiprows=1,
            usecols=(er_index,),
            dtype=str,
            unpack=True,
        )

        er_iter = np.array([er_data]).reshape(-1)
        er_iter = [int(s == "True") for s in er_iter]

        # except ValueError:
        #    print "Desired data not found in '%s'" %(sol_dir)
        #   continue

        plt.plot(e0_data, er_iter, ":ro")

    plt.xlabel("$E_0$")
    plt.ylabel("Ergo region")

    plt.show()

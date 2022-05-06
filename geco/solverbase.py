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
This module implements common functionality for the solvers.
"""

import math
import os
import sys
from os.path import join as pj

import numpy
from dolfin import *
#from dolfin import FunctionSpace, Parameters, Point
from mshr import *
from ufl.algorithms import extract_coefficients

from geco.anderson import *
from geco.models import *


def _dict2table(dict, title):
    "Convert dictionary to table"
    table = Table(title)
    keys = sorted(dict.keys())
    for key in keys:
        table.set(key, "value", dict[key])
    return table

class SolverBase:
    def __init__(self, solver_prefix):

        # Make FEniCS info print only on one processor
        #parameters.std_out_all_processes = False

        # Set up directories
        library_dir = os.path.dirname(os.path.abspath(__file__))
        geometry_dir = os.path.join(library_dir, "geometries")
        solution_dir = os.path.join("solutions", solver_prefix)

        # Create parameter set
        self.parameters = Parameters()

        # Discretization parameters
        self.parameters.add(Parameters("discretization"))
        self.parameters["discretization"].add("mass", 1.0)
        self.parameters["discretization"].add("domain_radius", 25)
        self.parameters["discretization"].add("maxiter", 1000)
        self.parameters["discretization"].add("theta", 1.0)
        self.parameters["discretization"].add("adaptive_theta", False)
        self.parameters["discretization"].add("anderson_depth", 3)
        self.parameters["discretization"].add("tolerance", 1e-3)
        self.parameters["discretization"].add("krylov_tolerance", 1e-9)
        self.parameters["discretization"].add("num_steps", 10)
        self.parameters["discretization"].add("degree", 1)
        self.parameters["discretization"].add("resolution", 32)
        self.parameters["discretization"].add("resolution_3d", 8)
        self.parameters["discretization"].add("mesh_prefix", "halfcircle")

        # Output parameters
        self.parameters.add(Parameters("output"))
        self.parameters["output"].add("library_directory", library_dir)
        self.parameters["output"].add("geometry_directory", geometry_dir)
        self.parameters["output"].add("solution_directory", solution_dir)
        self.parameters["output"].add("suffix", "")
        self.parameters["output"].add("plot_iteration", True)
        self.parameters["output"].add("plot_solution", False)
        self.parameters["output"].add("save_solution", True)
        self.parameters["output"].add("save_solution_3d", False)
        self.parameters["output"].add("save_point_cloud", False)
        self.parameters["output"].add("save_iterations", False)
        self.parameters["output"].add("save_residuals", False)

        # Override some parameters when --hires option is given
        if len(sys.argv) > 1 and sys.argv[1] == "--hires":
            info("Running simulation with high resolution")
            self.parameters["discretization"]["num_steps"] = 100
            self.parameters["discretization"]["num_steps"] = 50
            self.parameters["discretization"]["resolution"] = 512
            self.parameters["discretization"]["resolution_3d"] = 128
            self.parameters["discretization"]["tolerance"] = 1e-4
        if len(sys.argv) > 1:
            solution_dir = sys.argv[-1]
            self.parameters["output"]["solution_directory"] = solution_dir

        # Check whether plotting should be disabled
        self._dolfin_noplot = False
        if "DOLFIN_NOPLOT" in os.environ:
            info("Plotting disabled (DOLFIN_NOPLOT set).")

        # Initialize adaptive relaxation (using bisection)
        self._theta = 1.0
        self._theta_max = None
        self._theta_init = False

    def _init_solve(self):
        "Initializations at the beginning of solve call"

        # Get parameters
        solution_dir = self.parameters["output"]["solution_directory"]

        # Create file for saving density during iterations
        if self.parameters["output"]["save_iterations"]:
            self._density_file = XDMFFile(
                mpi_comm_world(), pj(solution_dir, "iterations.xdmf")
            )
            self._density_file.parameters.flush_output = True
        else:
            self._density_file = None

        # Create list of residuals
        self._residuals = []

    def _generate_mesh(self):

        # Note that we generate the mesh with unit radius and
        # then scale the mesh so as not to confuse the mshr
        # resolution...

        # Get parameters
        R = self.parameters["discretization"]["domain_radius"]
        N = self.parameters["discretization"]["resolution"]

        # Define domain (unit half disk)
        circle = Circle(Point(0, 0), 1)
        rectangle = Rectangle(Point(0, 0), Point(2, 2))
        domain = circle * rectangle

        # Generate uniform mesh
        mesh = generate_mesh(domain, N)

        # Scale mesh
        x = mesh.coordinates()
        x *= R

        return mesh

    def _generate_mesh_annulus(self, R, N):

        # Note that we generate the mesh with unit radius and
        # then scale the mesh so as not to confuse the mshr
        # resolution...

        # Get parameters
        R = self.parameters["discretization"]["domain_radius"]
        N = self.parameters["discretization"]["resolution"]

        # Define domain
        big_circle = Circle(Point(0, 0), 2)
        circle = Circle(Point(0, 0), 1)
        rectangle = Rectangle(Point(0, 0), Point(4, 4))
        domain = big_circle * rectangle - circle

        # Generate uniform mesh
        mesh = generate_mesh(domain, N)

        # Scale mesh
        x = mesh.coordinates()
        x *= R

        return mesh

    def _get_theta(self):
        "Try to find optimal theta"

        # If not adaptive theta, just return numerical value
        if not self.parameters["discretization"].adaptive_theta:
            return self.parameters["discretization"].theta

        # Check whether it's time to get started
        if not self._theta_init:

            # If we don't have enough residuals, keep going
            if len(self._residuals) < 5:
                return self._theta

            # If we haven't reached asymptotic regime, keep going
            r1, r2, r3 = self._residuals[-3:]
            if abs(r1 / r2 - r2 / r3) > 0.05:
                return self._theta

            info("Initializing adaptive relaxation")
            self._theta_init = True

        # If residual increases, set right end-point and decrease theta
        if self._residuals[-1] > self._residuals[-2]:
            self._theta_max = self._theta
            self._theta *= 0.5
            info("Residual increasing, setting theta_max = %.3g" % self._theta_max)

        # Find right end-point by carefully increasing theta
        elif self._theta_max is None:
            self._theta *= 1.05

        # Try approaching right end-point
        else:
            print(self._theta, self._theta_max)
            self._theta = (
                2.0 * self._theta * self._theta_max / (self._theta + self._theta_max)
            )

        info("theta = %.3g" % self._theta)

        return self._theta

    def _postprocess(
        self,
        ansatzes,
        solutions,
        flat_solutions,
        names,
        residual_functions,
        matter_components,
        matter_names,
    ):

        # File access sometimes fails in parallel and crashes the
        # solution, in particular the call to os.makedirs, so wrap
        # this in a try/except.

        try:
            self._print_data(ansatzes)
            self._save_solutions(solutions, names)
            self._save_solutions(matter_components, matter_names)
            self._save_residual_functions(residual_functions, names)
            self._save_flat(flat_solutions, names)
            self._save_solution_3d(solutions[-1])
            self._save_point_cloud(solutions[-1])
            self._plot_solutions(solutions[:-1], names[:-1])
            self._save_data()  # do this last as it may break
        except:
            warning("Postprocessing failed: %s" % str(sys.exc_info()[0]))

    def _print_data(self, ansatzes):
        "Pretty-print data"

        # Print ansatzes
        for ansatz in ansatzes:
            info("")
            info(ansatz.parameters, True)

        # Append all parameters (except output parameters)
        self.data.update(self.parameters["discretization"].to_dict())
        for ansatz in ansatzes:
            self.data.update(ansatz.parameters.to_dict())

        # Print data
        if MPI.rank(mpi_comm_world()) == 0:
            info("")
            info(_dict2table(self.data, "data"), True)
            info("")

    def _save_data(self):
        "Save data to file"

        # Do this only on processor 0
        if MPI.rank(mpi_comm_world()) > 0:
            return

        # Get parameters
        solution_dir = self.parameters["output"]["solution_directory"]

        # Create directory if it does not yet exist
        if not os.path.exists(solution_dir):
            os.makedirs(solution_dir)

        # Name of file
        filename = os.path.join(solution_dir, "data.csv")
        info("Appending data to file %s." % filename)

        # Get keys
        keys = sorted(self.data.keys())

        # Write header if not written before
        if not os.path.isfile(filename):
            f = open(filename, "w")
            f.write(",".join(str(k) for k in keys) + "\n")
            f.close()

        # Append data
        data_line = ""
        for key in keys:
            val = self.data[key]
            if isinstance(val, float):
                data_line += "%.16g" % val + ","
            else:
                data_line += str(val) + ","
        data_line = data_line.strip(",")

        f = open(filename, "a")
        f.write(data_line + "\n")
        f.close()

    def _save_solutions(self, solutions, names):
        "Save solutions to file"

        # Check whether to save solution
        if not self.parameters["output"]["save_solution"]:
            return

        # Get parameters
        R = self.parameters["discretization"]["domain_radius"]
        suffix = self.parameters["output"]["suffix"]
        solution_dir = self.parameters["output"]["solution_directory"]

        # Save mesh to XML format
        f = File(pj(solution_dir, "mesh.xml.gz"))
        f << solutions[0].function_space().mesh()

        # Save solutions to XML format
        for solution, name in zip(solutions, names):
            f = File(pj(solution_dir, "%s_%d%s.xml.gz" % (name, R, suffix)))
            f << solution

        # Save solutions to XDMF format
        for solution, name in zip(solutions, names):
            f = XDMFFile(
                mpi_comm_world(), pj(solution_dir, "%s_%d%s.xdmf" % (name, R, suffix))
            )
            f.write(solution)

    def _save_residual_functions(self, residual_functions, names):
        "Save residual functions to file"

        # Check whether to save residuals
        if not self.parameters["output"].save_residuals:
            return

        # Extract field names and solution directory
        field_names = [names[n] for n in xrange(4)]
        solution_dir = self.parameters["output"]["solution_directory"]

        # Save solutions to XDMF format
        for resf, fname in zip(residual_functions, field_names):
            f = XDMFFile(mpi_comm_world(), pj(solution_dir, "%s_residual.xdmf" % fname))
            f.write(resf)

    def _save_flat(self, solutions, names):
        "Save flat solutions on external annulus to file"

        # Check whether to save solution
        if not self.parameters["output"].save_solution:
            return

        # Get parameters
        R = self.parameters["discretization"]["domain_radius"]
        N = self.parameters["discretization"]["resolution"]
        suffix = self.parameters["output"]["suffix"]
        solution_dir = self.parameters["output"]["solution_directory"]

        # Generate solutions
        annulus_mesh = self._generate_mesh_annulus(R, N)
        Z = FunctionSpace(annulus_mesh, "Lagrange", 1)
        _solutions = [interpolate(solution, Z) for solution in solutions]

        # Save solutions to XDMF format
        for _solution, name in zip(_solutions, names):
            f = XDMFFile(
                mpi_comm_world(), pj(solution_dir, "%s_R_%d%s.xdmf" % (name, R, suffix))
            )
            f.write(_solution)

    def _save_solution_3d(self, RHO):
        "Save 3D solution to file"

        # Check whether to save solution
        if not self.parameters["output"].save_solution_3d:
            return

        # Get parameters
        R = self.parameters["discretization"]["domain_radius"]
        suffix = self.parameters["output"]["suffix"]
        solution_dir = self.parameters["output"]["solution_directory"]
        resolution_3d = self.parameters["discretization"].resolution_3d

        # Generate solution
        info("Computing 3D representation of density")
        n = resolution_3d
        box = BoxMesh(Point(-R, -R, -R), Point(R, R, R), n, n, n)
        V3D = FunctionSpace(box, "Lagrange", 1)
        rho3d = interpolate(Density3D(RHO), V3D)

    def _save_point_cloud(self, RHO):
        "Save point cloud to file"

        # Check whether to save solution
        if not self.parameters["output"].save_point_cloud:
            return

        # Get parameters
        R = self.parameters["discretization"]["domain_radius"]
        suffix = self.parameters["output"]["suffix"]
        solution_dir = self.parameters["output"]["solution_directory"]
        m = self.data["mass"]

        # Generate solution
        info("Computing point cloud representation of density")
        rho = PointCloud(RHO, R, m, 64, 50000)

        # Save solution
        filename = pj(solution_dir, "point_cloud_%d%s.xdmf" % (R, suffix))
        rho.save_data(filename)

    def _save_density(self, RHO, iter):
        "Save density during iterations"

        # Check whether to save density
        if not self.parameters["output"]["save_iterations"]:
            return

        # Save density to file
        self._density_file.write(RHO, float(iter))

    def _save_residual(self, residual):
        "Save residual to file"

        # Store residual
        self._residuals.append(residual)

        # Get parameters
        solution_dir = self.parameters["output"]["solution_directory"]

        # Create directory if it does not yet exist
        if not os.path.exists(solution_dir):
            os.makedirs(solution_dir)

        # Write residuals
        filename = pj(solution_dir, "residuals.csv")
        with open(filename, "w") as f:
            f.write(",".join("%.16g" % r for r in self._residuals))

    def _plot_solutions(self, solutions, names):
        "Plot solution"

        # Check whether to plot solution
        if self._dolfin_noplot or not self.parameters["output"]["plot_solution"]:
            return

        # Plot solutions
        for solution, name in zip(solutions, names):
            plot(solution, title=name)
        interactive()

    def _plot_density(self, RHO):
        "Plot density during iterations"

        # Check whether to plot iteration
        if self._dolfin_noplot or not self.parameters["output"]["plot_solution"]:
            return

        # Plot density
        plot(RHO, title="Density")

"This module implements common functionality for the solvers."

import os, sys, math, numpy
from os.path import join as pj

from dolfin import *
from ufl.algorithms import extract_coefficients

from models import *

def _dict2table(dict, title):
    "Convert dictionary to table"
    table = Table(title)
    keys = sorted(dict.keys())
    for key in keys:
        table.set(key, "value", dict[key])
    return table

class SolverBase:

    def __init__(self):

        # Set up directories
        library_dir  = os.path.dirname(os.path.abspath(__file__))
        geometry_dir = os.path.join(library_dir, "geometries")
        solution_dir = os.path.join("solutions", "ev")

        # Create parameter set
        self.parameters = \
            Parameters(discretization = Parameters("discretization"),
                       output = Parameters("output"))

        # Discretization parameters
        self.parameters.discretization.add("mass", 1.0)
        self.parameters.discretization.add("radius", 25)
        self.parameters.discretization.add("maxiter", 1000)
        self.parameters.discretization.add("theta", 1.0)
        self.parameters.discretization.add("tolerance", 1e-3)
        self.parameters.discretization.add("krylov_tolerance", 1e-9)
        self.parameters.discretization.add("num_steps", 10)
        self.parameters.discretization.add("resolution", 32)
        self.parameters.discretization.add("resolution_3d", 8)
        self.parameters.discretization.add("mesh_prefix", "halfcircle")

        # Output parameters
        self.parameters.output.add("library_directory",  library_dir)
        self.parameters.output.add("geometry_directory", geometry_dir)
        self.parameters.output.add("solution_directory", solution_dir)
        self.parameters.output.add("suffix",             "")
        self.parameters.output.add("plot_iteration",     True)
        self.parameters.output.add("plot_solution",      False)
        self.parameters.output.add("save_solution",      True)
        self.parameters.output.add("save_solution_3d",   False)
        self.parameters.output.add("save_point_cloud",   False)

        # Override some parameters when --hires option is given
        if len(sys.argv) > 1 and sys.argv[1] == "--hires":
            info("Running simulation with high resolution")
            self.parameters.discretization.num_steps = 100
            self.parameters.discretization.num_steps = 50
            self.parameters.discretization.resolution = 512
            self.parameters.discretization.resolution_3d = 128
            self.parameters.discretization.tolerance = 1e-4
        if len(sys.argv) > 1:
            solution_dir = sys.argv[-1]
            self.parameters.output.solution_directory = solution_dir

    def _read_mesh(self):

        # Get parameters
        R = self.parameters.discretization.radius
        N = self.parameters.discretization.resolution
        geometry_dir = self.parameters.output.geometry_directory
        mesh_prefix = self.parameters.discretization.mesh_prefix

        # Read from file
        mesh_file = pj(geometry_dir, "%s_%d_%d_mesh.xml.gz" % (mesh_prefix, R, N))
        mesh = Mesh(mesh_file)

        return mesh

    def _postprocess(self, ansatzes, solutions, flat_solutions, names):

        self._print_data(ansatzes)
        self._save_data()
        self._save_solutions(solutions, names)
        self._save_flat(flat_solutions, names)
        self._save_solution_3d(solutions[-1])
        self._save_point_cloud(solutions[-1])
        self._plot_solutions(solutions[:-1], names[:-1])

    def _print_data(self, ansatzes):
        "Pretty-print data"

        # Print ansatzes
        for ansatz in ansatzes:
            info("")
            info(ansatz.parameters, True)

        # Append all parameters (except output parameters)
        self.data.update(self.parameters.discretization.to_dict())
        for ansatz in ansatzes:
            self.data.update(ansatz.parameters.to_dict())

        # Print data
        info("")
        info(_dict2table(self.data, "data"), True)
        info("")

    def _save_data(self):
        "Save data to file"

        # Get parameters
        solution_dir = self.parameters.output.solution_directory

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
        data_line=""
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
        if not self.parameters.output.save_solution: return

        # Get parameters
        R = self.parameters.discretization.radius
        suffix = self.parameters.output.suffix
        solution_dir = self.parameters.output.solution_directory

        # Save solutions to XML format
        for solution, name in zip(solutions, names):
            f = File(pj(solution_dir, "%s_%d%s.xml" % (name, R, suffix)))
            f << solution

        # Save solutions in VTK format
        for solution, name in zip(solutions, names):
            f = File(pj(solution_dir, "%s_%d%s.pvd" % (name, R, suffix)))
            f << solution

    def _save_flat(self, solutions, names):
        "Save flat solutions on external annulus to file"

        # Check whether to save solution
        if not self.parameters.output.save_solution: return

        # Get parameters
        R = self.parameters.discretization.radius
        N = self.parameters.discretization.resolution
        suffix = self.parameters.output.suffix
        geometry_dir = self.parameters.output.geometry_directory
        solution_dir = self.parameters.output.solution_directory

        # Generate solutions
        f = pj(geometry_dir, "halfannulus_%d_%d_mesh.xml.gz" % (R, N))
        annulus_mesh = Mesh(f)
        Z = FunctionSpace(annulus_mesh, "Lagrange", 1)
        _solutions = [interpolate(solution, Z) for solution in solutions]

        # Save solutions
        for _solution, name in zip(_solutions, names):
            f = File(pj(solution_dir, "%s_R_%d%s.pvd" % (name, R, suffix)))
            f << _solution

    def _save_solution_3d(self, RHO):
        "Save 3D solution to file"

        # Check whether to save solution
        if not self.parameters.output.save_solution_3d: return

        # Get parameters
        R = self.parameters.discretization.radius
        suffix = self.parameters.output.suffix
        solution_dir = self.parameters.output.solution_directory
        resolution_3d = self.parameters.discretization.resolution_3d

        # Generate solution
        info("Computing 3D representation of density")
        n = resolution_3d
        box = BoxMesh(Point(-R, -R, -R), Point(R, R, R), n, n, n)
        V3D = FunctionSpace(box, "Lagrange", 1)
        rho3d = interpolate(Density3D(RHO), V3D)

    def _save_point_cloud(self, RHO):
        "Save point cloud to file"

        # Check whether to save solution
        if not self.parameters.output.save_point_cloud: return

        # Get parameters
        R = self.parameters.discretization.radius
        suffix = self.parameters.output.suffix
        solution_dir = self.parameters.output.solution_directory
        m = self.data["mass"]

        # Generate solution
        info("Computing point cloud representation of density")
        rho = PointCloud(RHO, R, m, 64, 50000)

        # Save solution
        filename = pj(solution_dir, "point_cloud_%d%s.xdmf" % (R, suffix))
        rho.save_data(filename)

    def _plot_solutions(self, solutions, names):
        "Plot solution"

        # Check whether to plot solution
        if not self.parameters.output.plot_solution: return

        # Plot solutions
        for solution, name in zip(solutions, names):
            plot(solution, title=name)
        interactive()
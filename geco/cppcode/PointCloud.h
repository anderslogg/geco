// C++ code for point cloud representation of density. Note that this
// is really not an Expression and not used as such but Expression
// subclassing is used as a convenient way to inject C++ code in
// Python.

// ---------------------------------------------------------
// Copyright 2019 Anders Logg, Ellery Ames, Håkan Andréasson

// This file is part of GECo.
// GECo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// GECo is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with GECo. If not, see <https://www.gnu.org/licenses/>.

#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/geometry/Point.h>
#include "dolfin/io/XDMFFile.h"
#include <random>

class PointCloud : public dolfin::Expression
{
public:
  // Constructor
  PointCloud() : Expression(),
                 _R(0), _resolution(0), _num_points(0) {}

  // Evaluation
  void eval(Eigen::Ref<Eigen::VectorXd> values, 
            Eigen::Ref<const Eigen::VectorXd> x,
            const ufc::cell &cell) const
  {
    dolfin::error("Evaluation of point cloud not implemented.");
  }

  // Set parameters
  void set_parameters(std::shared_ptr<const dolfin::Function> rho,
                      double R,
                      double M,
                      std::size_t resolution,
                      std::size_t num_points)
  {
    _rho = rho;
    _R = R;
    _M = M;
    _resolution = resolution;
    _num_points = num_points;
  }

  // Save point cloud to file
  void save_data(std::string filename)
  {
    dolfin_assert(_rho);
    dolfin_assert(_R > 0.0);
    dolfin_assert(_M > 0.0);
    dolfin_assert(_resolution > 0);
    dolfin_assert(_num_points > 0);

#ifdef HAS_HDF5

    // Extract parameters
    const std::size_t n = _resolution;
    const std::size_t N = _num_points;
    std::size_t empty_boxes = 0;

    // Compute some data for point cloud
    const double h = 2.0 * _R / static_cast<double>(n);
    dolfin::info("Number of points:         %d", N);
    dolfin::info("Size of sample area:      %g", 2. * _R);
    dolfin::info("Number of sample boxes:   %d", n * n * n);
    dolfin::info("Size of sample boxes:     %g", h);

    // Create empty point cloud data
    std::vector<dolfin::Point> point_cloud;
    std::vector<double> point_cloud_values;

    // Initialize distribution used to place points in sample boxes
    std::random_device random_seed_factory;
    std::mt19937 random_generator(random_seed_factory());
    std::uniform_real_distribution<double> uniform_point_distribution(0., 1.);

    // Iterate over boxes
    for (std::size_t nx = 0; nx < n; nx++)
    {
      const double x0 = -_R + nx * h;
      for (std::size_t ny = 0; ny < n; ny++)
      {
        const double y0 = -_R + ny * h;
        for (std::size_t nz = 0; nz < n; nz++)
        {
          const double z0 = -_R + nz * h;

          // Note: We flip the axes here to get a plot that
          // can be viewed at the same time as the xy-plane plot
          // of the density in Paraview.

          // Note: For some strange reason abs() does not return
          // anything sensible (double/int issue?) so we use our own.

          // Compute cylindrical coordinates for midpoint of box
          const double _x = x0 + 0.5 * h;
          const double _y = y0 + 0.5 * h;
          const double _z = z0 + 0.5 * h;
          const double s = sqrt(_x * _x + _z * _z);
          const double z = (_y >= 0.0 ? _y : -_y);

          // Skip if point is outside of domain (sphere)
          if (sqrt(_x * _x + _y * _y + _z * _z) > 0.999 * _R)
          {
            empty_boxes += 1;
            continue;
          }

          // Evaluate density at midpoint of box
          const double box_density = (*_rho)(s, z);

          // Compute number of points in box
          const double box_mass = box_density * h * h * h;
          std::size_t points_in_box = static_cast<std::size_t>(box_mass * N / _M);

          // Randomize positions of points in box
          for (std::size_t i = 0; i < points_in_box; i++)
          {
            // Generate a random point in the box
            const double rx = uniform_point_distribution(random_generator);
            const double ry = uniform_point_distribution(random_generator);
            const double rz = uniform_point_distribution(random_generator);
            const double x = x0 + rx * h;
            const double y = y0 + ry * h;
            const double z = z0 + rz * h;
            
            // Add box to point cloud
            dolfin::Point p(x, y, z);
            point_cloud.push_back(p);
            point_cloud_values.push_back(box_density);
          }
        }
      }
    }

    // Save to file
    dolfin::info("Number of sample boxes populated %d", n*n*n - empty_boxes);
    dolfin::info("Total number of points:          %d", point_cloud.size());
    dolfin::XDMFFile file(_rho->function_space()->mesh()->mpi_comm(), filename);
    file.write(point_cloud, point_cloud_values);

#elif

    warning("Unable to save point cloud data: missing HDF5 library.");

#endif
  }

private:
  // Axially symmetric density (2D)
  std::shared_ptr<const dolfin::Function> _rho;

  // Parameters
  double _R; // size of box [-R, R]^3
  double _M; // total mass (integral of RHO)
  std::size_t _resolution;
  std::size_t _num_points;
};
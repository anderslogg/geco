// C++ code for 3D representation of density

// ---------------------------------------------------------
// Copyright 2019 Anders Logg, Ellery Ames, Håkan Andréasson

// This file is part of GECo.
// GECo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// GECo is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with GECo. If not, see <https://www.gnu.org/licenses/>.

#ifndef DENSITY3D_H
#define DENSITY3D_H

#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <iostream>

class Density3D : public dolfin::Expression
{
public:
  // Constructor
  Density3D() : Expression() {}

  // Evaluation
  void eval(Eigen::Ref<Eigen::VectorXd> values, 
            Eigen::Ref<const Eigen::VectorXd> x,
            const ufc::cell &cell) const
  {
    dolfin_assert(_rho);

    // Note: We flip the axes here to get a plot that
    // can be viewed at the same time as the xy-plane plot
    // of the density in Paraview.

    // Note: For some strange reason abs() does not return
    // anything sensible (double/int issue?) so we use our own.

    // Note: Seem to need a small eps here even if allow_extrapolation
    // is turned on. Strange, but doesn't matter.

    // Extract cylindrical coordinates
    const double eps = 1e-6;
    const double s = sqrt(x[0] * x[0] + x[2] * x[2]) + eps;
    const double z = (x[1] >= 0.0 ? x[1] : -x[1]);

    // std::cout << "s = " << s << " z = " << z << std::endl;

    // Evaluate at point
    values[0] = (*_rho)(s, z);
  }

  // Set density
  void set_density(std::shared_ptr<const dolfin::Function> rho)
  {
    _rho = rho;
  }

private:
  // Axially symmetric density (2D)
  std::shared_ptr<const dolfin::Function> _rho;
};

#endif
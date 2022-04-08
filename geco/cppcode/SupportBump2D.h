// C++ code for 2D bump function on support of energy density

// ---------------------------------------------------------
// Copyright 2019 Anders Logg, Ellery Ames, Håkan Andréasson

// This file is part of GECo.
// GECo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

// GECo is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with GECo. If not, see <https://www.gnu.org/licenses/>.

class SupportBump2D : public Expression
{
public:
  // Constructor
  SupportBump2D() : Expression() {}

  // Evaluation
  void eval(Array<double> &values, const Array<double> &x,
            const ufc::cell &cell) const
  {
    dolfin_assert(_rho);

    // Extract cylindrical coordinates
    const double s = x[0];
    const double z = x[1];

    // Evaluate at point
    const double density = (*_function)(s, z);

    if (density > 0.0001)
      values[0] = 1.0;
    else
      values[0] = 0.0;
  }

  // Set density
  void set_density(std::shared_ptr<const Function> rho)
  {
    _rho = rho;
  }

private:
  // Axially symmetric function (2D)
  std::shared_ptr<const Function> _function;
};

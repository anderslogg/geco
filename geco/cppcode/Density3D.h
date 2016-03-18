// C++ code for 3D representation of density

class Density3D : public Expression
{
public:

  // Constructor
  Density3D() : Expression() {}

  // Evaluation
  void eval(Array<double>& values, const Array<double>& x,
            const ufc::cell& cell) const
  {
    dolfin_assert(_rho);

    // Note that we flip the axes here to get a plot that
    // can be viewed at the same time as the xy-plane plot
    // of the density in Paraview.

    // Extract cylindrical coordinates
    const double s = sqrt(x[0]*x[0] + x[2]*x[2]);
    const double z = x[1];

    // Evaluate at point
    values[0] = (*_rho)(s, z);
  }

  // Set density
  void set_density(std::shared_ptr<const Function> rho)
  {
    _rho = rho;
  }

private:

  // Axially symmetric density (2D)
  std::shared_ptr<const Function> _rho;

};

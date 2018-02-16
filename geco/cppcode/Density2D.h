// C++ code for 2D representation of density

class Density2D : public Expression
{
public:

  // Constructor
  Density2D() : Expression() {}

  // Evaluation
  void eval(Array<double>& values, const Array<double>& x,
            const ufc::cell& cell) const
  {
    dolfin_assert(_rho);

    // Note: Seem to need DOLFIN_EPS here even if allow_extrapolation
    // is turned on. Strange, but doesn't matter.

    // Extract cylindrical coordinates
    const double s = (x[0] >= 0.0 ? x[0] : -x[0]) + DOLFIN_EPS;
    const double z = (x[1] >= 0.0 ? x[1] : -x[1]);

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

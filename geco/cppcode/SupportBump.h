// C++ code for bump function on support of energy density

class SupportBump : public Expression
{
public:

  // Constructor
  SupportBump() : Expression() {}

  // Evaluation
  void eval(Array<double>& values, const Array<double>& x,
            const ufc::cell& cell) const
  {
    dolfin_assert(_rho);

    // Extract cylindrical coordinates
    const double s = x[0];
    const double z = x[1];

    // Evaluate at point
    const double density = (*_rho)(s, z);

    if (density > 0.0001)
      values[0] = 1;
    else
      values[0] = 0;
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

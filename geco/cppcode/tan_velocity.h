// C++ template code for density functionals. This defines a (sort of)
// base class, with functions to be defined by specific ansatzes. This
// is used instead of C++ inheritance since that is not supported by
// the DOLFIN JIT compiler.

namespace dolfin
{

  class VPAnsatz : public Expression
  {
  public:

    // Constructor
    VPAnsatz() : Expression(), _resolution(0)
    {
      // Set default parameter values
      init_parameters();
    }

    // Note: We need to implement both eval functions (with and
    // without the cell argument). The first one is used during
    // assembly (for speed), while the second will be used to compute
    // a 3D representation. For the 3D representation, the cell
    // argument is unknown so that step will involve building a search
    // tree.

    // Evaluate at given point in cell
    void eval(Array<double>& values,
              const Array<double>& x,
              const ufc::cell& cell) const
    {
      // Evaluate potential at point in cell
      dolfin_assert(_U);
      _U->eval(values, x, cell);
      const double U = values[0];

      // Evaluate ansatz
      values[0] = eval(U, x);
    }

    // Evaluate at given point
    void eval(Array<double>& values,
              const Array<double>& x) const
    {
      // Evaluate potential at point
      dolfin_assert(_U);
      _U->eval(values, x);
      const double U = values[0];

      // Evaluate ansatz
      values[0] = eval(U, x);
    }

    // Evaluation of ansatz
    // FIXME: Evaluate each species ansatz seperately
    double eval(double U,
                const Array<double>& x) const
    {
      // Get coordinates
      const double rho = x[0];
      const double z = x[1];

      // Check cut-off
      if (U >= E0)
        return 0.0;

      // Get resolution
      const double n = _resolution;
      dolfin_assert(n > 0);

      // Reset integral
      double I = 0.0;

      // Compute integration limits for E-integral
      const double Ea = U;
      const double Eb = E0;
      const double dE = (Eb - Ea) / static_cast<double>(n);
      if (dE <= 0.0) error("Strange stuff: dE <= 0.0");

      // Integrate over E
      for (std::size_t i = 0; i < n; i++)
      {
        // Compute integration variable
        const double E = Ea + static_cast<double>(i)*dE + 0.5*dE;

        // Check expression for integration limit
        if (E < U) error("Strange stuff: E < U");

        // Compute integration limits for s-integral
        const double sm = std::sqrt(2.0*(E - U));
        const double sa = -sm;
        const double sb = sm;

        // Compute step size for s-integral
        const double ds = (sb - sa) / static_cast<double>(n);
        if (ds < 0.0) error("Strange stuff: ds < 0.0");

        // Integrate over s
        for (std::size_t j = 0; j < n; j++)
        {
          // Compute integration variable
          const double s = sa + static_cast<double>(j)*ds + 0.5*ds;

          // Evaluate ansatz and add to integral
          const double L = rho*s;

          I += s*ansatz(E, L)*ds*dE;
        }
      }

      // Scale integral
      I *= 2.0*DOLFIN_PI;

      // Update radius of support
      const double R = std::sqrt(x[0]*x[0] + x[1]*x[1]);
      if (I > DOLFIN_EPS && R > _R)
        _R = R;

      return I;
    }

    // Set potential
    void set_fields(std::shared_ptr<const Function> U)
    {
      _U = U;
    }

    // Set integration parameters
    void set_integration_parameters(std::size_t resolution)
    {
      _resolution = resolution;
    }

    // Reset computation of radius of support
    void reset()
    {
      _R = 0;
    }

    // Return radius of support
    double radius_of_support()
    {
      return _R;
    }

    // The parameters
    Parameters parameters;

    // Member functions (to be defined by specific ansatz)
    %(member_functions)s

  private:

    // Number of steps to use in numerical integration
    std::size_t _resolution;

    // Radius of support
    mutable double _R;

    // The potential
    std::shared_ptr<const Function> _U;

    // Member variables (to be defined by specific ansatz)
    %(member_variables)s

  };

}

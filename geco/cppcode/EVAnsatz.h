// C++ template code for density functionals. This defines a (sort of)
// base class, with functions to be defined by specific ansatzes. This
// is used instead of C++ inheritance since that is not supported by
// the DOLFIN JIT compiler.

namespace dolfin
{

  class EVAnsatz : public Expression
  {
  public:

    // Constructor
    EVAnsatz() : Expression(5), _resolution(0)
    {
      // Set default parameter values
      init_parameters();

      // Add rotation parameter
      parameters.add("rotation", true);

      // Add particle mass parameter
      parameters.add("particle_mass", 1.0);
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
      dolfin_assert(_NU);
      dolfin_assert(_BB);
      dolfin_assert(_MU);
      dolfin_assert(_WW);

      // Create value array for evaluation of fields. Note that we
      // cannot reuse the values array since it has length 3 and
      // the fields we evaluate are scalar.
      Array<double> _values(1);

      // Evaluate NU
      _NU->eval(_values, x, cell);
      const double NU = _values[0];

      // Evaluate BB
      _BB->eval(_values, x, cell);
      const double BB = _values[0];

      // Evaluate MU
      _MU->eval(_values, x, cell);
      const double MU = _values[0];

      // Evaluate WW
      _WW->eval(_values, x, cell);
      const double WW = _values[0];

      // Evaluate ansatz
      eval(values, NU, BB, MU, WW, x);
    }

    // Evaluate at given point
    void eval(Array<double>& values,
              const Array<double>& x) const
    {
      dolfin_assert(_NU);
      dolfin_assert(_BB);
      dolfin_assert(_MU);
      dolfin_assert(_WW);

      // Create value array for evaluation of fields. Note that we
      // cannot reuse the values array since it has length 3 and
      // the fields we evaluate are scalar.
      Array<double> _values(1);

      // Evaluate NU
      _NU->eval(_values, x);
      const double NU = _values[0];

      // Evaluate BB
      _BB->eval(_values, x);
      const double BB = _values[0];

      // Evaluate MU
      _MU->eval(_values, x);
      const double MU = _values[0];

      // Evaluate WW
      _WW->eval(_values, x);
      const double WW = _values[0];

      // Evaluate ansatz
      eval(values, NU, BB, MU, WW, x);
    }

    // Evaluation of ansatz
    void eval(Array<double>& values,
	      double NU, double BB, double MU, double WW,
	      const Array<double>& x) const
    {
      // We want WW non-negative
      _min_WW = std::min(_min_WW, WW);
      WW = std::max(0.0, WW);

      // Get coordinates
      const double rho = x[0];
      const double z = x[1];

      // Compute expressions independent of integration variables
      const double pi = DOLFIN_PI;
      const double K1 = exp(NU);

      // Check cut-off
      if (_particle_mass*K1 >= E0)
      {
        values[0] = 0.0;
        values[1] = 0.0;
        values[2] = 0.0;
        values[3] = 0.0;
        values[4] = 0.0;
        return;
      }

      // Get resolution
      const double n = _resolution;
      dolfin_assert(n > 0);

      // Reset integrals
      double I0 = 0.0;
      double I1 = 0.0;
      double I2 = 0.0;
      double I3 = 0.0;
      double I4 = 0.0;

      // Compute integration limits for h-integral
      const double ha = _particle_mass*K1;
      const double hb = E0;
      const double dh = (hb - ha) / static_cast<double>(n);
      if (dh <= 0.0) error("Strange stuff: dh <= 0.0");

      // Integrate over h
      for (std::size_t i = 0; i < n; i++)
      {
        // Compute integration variable
        const double h = ha + static_cast<double>(i)*dh + 0.5*dh;

        // Check expression for integration limit
        const double K2 = h*h / (K1*K1) - _particle_mass*_particle_mass;
        if (K2 < 0.0) error("Strange stuff: K2 < 0.0");

        // Compute integration limits for s-integral
        double sa = 0.0;
        const double sb = (BB/K1)*sqrt(K2);

        // Change lower limit if not rotating
        if (!_rotation)
          sa = -sb;

        // Compute step size for s-integral
        const double ds = (sb - sa) / static_cast<double>(n);
        if (ds < -10.0*DOLFIN_EPS)
	  {
	    std::cout << "BB:" << BB << "\n";
	    std::cout << "sb:" << sb << "\n";
	    std::cout << "ds:" << ds << "\n";
	    error("Strange stuff: ds < 0.0");
	  }

        // Integrate over s
        for (std::size_t j = 0; j < n; j++)
        {
          // Compute integration variable
          const double s = sa + static_cast<double>(j)*ds + 0.5*ds;

          // Compute variables for ansatz
          const double L = rho*s;
          const double E = h + WW*L;
          if (L < 0.0 && _rotation) error("Strange stuff: L < 0.0");
          if (E < 0.0) error("Strange stuff: E < 0.0");

          // Evaluate ansatz
          const double f = ansatz(E, L);
          if (f < 0.0) error("Strange stuff: f < 0.0");

          // Compute integrals
          I0 += E*E*f*ds*dh;
          I1 += (sb*sb - s*s)*f*ds*dh;
          I2 += s*s*f*ds*dh;
          I3 += s*E*f*ds*dh;
          I4 += f*h*ds*dh;
        }
      }

      // Compute XI
      const double XI = MU + NU;

      // Scale integrals to compute Phi00, Phi11, Phi33, Phi03, density distribution
      values[0] =  2.0*pi*exp(2.0*XI)*std::pow(BB, -1)*std::pow(K1, -4.0)*I0;     // Integral for PHI_00
      values[1] =  2.0*pi*exp(2.0*XI)*std::pow(BB, -3)*I1;                        // Integral for PHI_11
      values[2] =  2.0*pi*exp(2.0*XI)*std::pow(BB, -3)*I2;                        // Integral for PHI_33
      values[3] = -2.0*pi*exp(2.0*XI)*std::pow(BB, -1)*rho*I3;                    // Integral for PHI_03
      values[4] =  2.0*pi*exp(2.0*MU - 2.0*NU)*_particle_mass*I4;                 // Integral for the rest density

      // Update radius of support
      const double R = std::sqrt(x[0]*x[0] + x[1]*x[1]);
      if (I4 > DOLFIN_EPS && R > _R)
        _R = R;
    }

    // Set fields
    void set_fields(std::shared_ptr<const Function> NU,
		    std::shared_ptr<const Function> BB,
                    std::shared_ptr<const Function> MU,
		    std::shared_ptr<const Function> WW)
    {
      _NU = NU;
      _BB = BB;
      _MU = MU;
      _WW = WW;
    }

    // Reset computation of radius of support
    void reset()
    {
      _R = 0;
      _min_WW = 0;
    }

    // Return radius of support
    double radius_of_support() const
    {
      return _R;
    }

    // Return minimum value of WW
    double min_WW() const
    {
      return _min_WW;
    }

    // Set integration parameters
    void set_integration_parameters(std::size_t resolution)
    {
      _resolution = resolution;
      _rotation = parameters["rotation"];
      _particle_mass = parameters["particle_mass"];
    }

    // The parameters
    Parameters parameters;

    // Member functions (to be defined by specific ansatz)
    %(member_functions)s

  private:

    // Number of steps to use in numerical integration
    std::size_t _resolution;

    // Rotating or not
    bool _rotation;

    // Particle Mass
    double _particle_mass;

    // Radius of support
    mutable double _R;

    // Minimum value of WW (should not be negative)
    mutable double _min_WW;

    // The fields
    std::shared_ptr<const Function> _NU;
    std::shared_ptr<const Function> _BB;
    std::shared_ptr<const Function> _MU;
    std::shared_ptr<const Function> _WW;

    // Member variables (to be defined by specific ansatz)
    %(member_variables)s

  };

}

// Ansatz for Vlasov-Poisson with Evans type ansatz in E and L.

#include "VPAnsatz.h"


class VPEvansLPoly : public VPAnsatz
{
 public:
  VPEvansLPoly() : VPAnsatz()
  {
    // Set default parameter values
    init_parameters();
  }

  // Member functions

  double ansatz(double E, double L) const
  {
    if (E >= E0)
      return 0.0;

    //  if (std::abs(L) <= L0)
      //    return 0.0;

    return c*std::exp(E*std::pow(s0, -2.0))*std::pow(abs(L), 2.0) + std::exp(E*std::pow(v0, -2.0));
  }

  void init_parameters()
  {
    parameters.add("E0", -0.1);
    parameters.add("v0",  0.0);
    parameters.add("s0",  1.0);
    parameters.add("c",   0.1);

  }

  void read_parameters()
  {
    E0 = parameters["E0"];
    v0 = parameters["v0"];
    s0 = parameters["s0"];
    c  = parameters["c"];
  }

private:
  // Member variables
  // Note that the energy cutoff E0 is declared within VPAnsatz.

  double v0;
  double s0;
  double c;

}; // end class VPEvansLPoly
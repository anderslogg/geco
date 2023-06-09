// Polytropic ansatz for Vlasov-Poisson

#include "VPAnsatz.h"

class VPEPolyLPoly : public VPAnsatz
{
public:
  // Constructor calls base class constructor
  VPEPolyLPoly() : VPAnsatz() 
  {
      // Set default parameter values
      init_parameters();
  };

  // Member functions 
  
  double ansatz(double E, double L) const override
  {
    if (E0 <= E)
      return 0.0;

    if (std::abs(L) < L0)
      return 0.0;

    return std::pow(E0 - E, k)*std::pow(std::abs(L) - L0, l);
  }

  void init_parameters()
  {
    parameters.add("E0", -0.1);
    parameters.add("L0",  0.0);
    parameters.add("k",   0.0);
    parameters.add("l",   0.0);
  }

  void read_parameters()
  {
    E0 = parameters["E0"];
    L0 = parameters["L0"];
    k  = parameters["k"];
    l  = parameters["l"];
  }

private:
  // Member variables
  // Note that the energy cutoff E0 is declared within VPAnsatz.

  double L0;
  double k;
  double l;
}; // end class VPEPolyLPoly

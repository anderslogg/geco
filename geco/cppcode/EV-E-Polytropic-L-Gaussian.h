//  EV-E-Polytropic-L-Gaussian
//  --------------------------
//    Class with ansatz function
//    
//      f(E, L) = (E0 - E)**k * (1/L0) * exp( +/- (L/L0)**2 )
//    for 
//      E > E0    &&   L > 0 if rotating.
//  --------------------------

#include "EVAnsatz.h"

class EVEPolyLGauss : public EVAnsatz
{
public: 
  EVEPolyLGauss() : EVAnsatz()
  {
    init_parameters();
  }

  double ansatz(double E, double L) const
  {
    if (E0 <= E)
      return 0.0;

    if (_rotation && L < 0.0)
      return 0.0;

    return std::pow(E0 - E, k)*(1.0 / L0)*std::exp(pm*std::pow(L / L0, 2.0));
  }

  void init_parameters()
  {
    parameters.add("E0", 0.9);
    parameters.add("k",  0.0);
    parameters.add("L0", 0.1);
    parameters.add("pm", 1.0);
  }

  void read_parameters()
  {
    E0  = parameters["E0"];
    k   = parameters["k"];
    L0  = parameters["L0"];
    pm  = parameters["pm"];
  } 

private:
  // Member variables

  double k;
  double L0;
  double pm;
}; // end of class EVEPolyLGauss
// Member functions


#include "EVAnsatz.h"

class EVEPolyLPoly : public EVAnsatz
{
public:
  EVEPolyLPoly() : EVAnsatz()
  {
    init_parameters();
  }

  double ansatz(double E, double L) const
  {
    if (E0 <= E)
      return 0.0;

    if (!_rotation)
      L = std::abs(L);

    if (L < L0)
      return 0.0;

    return std::pow(E0 - E, k)*std::pow(L - L0, l);
  }

  void init_parameters()
  {
    parameters.add("E0", 0.9);
    parameters.add("L0", 0.0);
    parameters.add("k",  0.0);
    parameters.add("l",  0.0);
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

  //double E0;
  double L0;
  double k;
  double l;
}; // end class EVEPolyLPoly
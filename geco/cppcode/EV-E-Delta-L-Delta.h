// Member functions

#include "EVAnsatz.h"

class EVEDeltaLDelta : public EVAnsatz
{
public:
  EVEDeltaLDelta() : EVAnsatz()
  {
    init_parameters();
  }

  double ansatz(double E, double L) const
  {
    if (E >= E0)
      return 0.0;

    if (E <= E0-e)
      return 0.0;

    if (!_rotation)
      L = std::abs(L);

    if (L <= L0 - l)
      return 0.0; 

    if (L >= L0 + l)
      return 0.0; 

    return 1.0/(4*e*l);
  }

  void init_parameters()
  {
    parameters.add("E0", 0.9);
    parameters.add("L0", 0.0);
    parameters.add("e",  0.05);
    parameters.add("l",  0.05);
  }

  void read_parameters()
  {
    E0 = parameters["E0"];
    L0 = parameters["L0"];
    e  = parameters["e"];
    l  = parameters["l"];
  }

private:
  // Member variables

  double L0;
  double e;
  double l;
}; // end class EVEDeltaLDelta
// Member functions

#include "EVAnsatz.h"

class EVEPolyLAndreasson : public EVAnsatz
{
public:
  EVEPolyLAndreasson() : EVAnsatz()
  {
    init_parameters();
  }

  double ansatz(double E, double L) const
  {
    if (E0 <= E)
      return 0.0;

    if (_rotation && L <= 0.0)
      return 0.0;

    if (!_rotation)
      L = std::abs(L);

    if (L >= 1.0 / Q)
      return 0.0;

    return std::pow(E0 - E, k)*std::pow(1 - Q*L, l);
  }

  void init_parameters()
  {
    parameters.add("E0", 0.9);
    parameters.add("k",  0.0);
    parameters.add("l",  0.0);
    parameters.add("Q",  0.0);
  }

  void read_parameters()
  {
    E0 = parameters["E0"];
    k  = parameters["k"];
    l  = parameters["l"];
    Q  = parameters["Q"];
  }

private:
  // Member variables

  double k;
  double l;
  double Q;
}; // end class EVEPolyLAndreasson
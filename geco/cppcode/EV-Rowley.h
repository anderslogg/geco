// Member functions

double chi(double E, double L) const
{
  // FIXME: Must look up this ansatz if used and check the following expression.
  return c0 - E + W*L - 0.5*std::pow(L / r0, 2.0);
}

double ansatz(double E, double L) const
{
  const double _chi = chi(E, L);
  
  if ( _chi <= 0)
    return 0.0;

  /* if (E >= E0) */
  /*   return 0.0; */

  if (_rotation && L <= 0.0)
    return 0.0;

  return std::exp(b*_chi);
  //  return std::pow(_chi, b);
  //  return std::log(1.0 + _chi);
  //  return 1.0/(1.0 + _chi);
}

void init_parameters()
{
  parameters.add("E0", 0.9);
  parameters.add("c0", 0.9);
  parameters.add("W",  1.0);
  parameters.add("r0", 1.0);
  parameters.add("b",  0.1);
}

void read_parameters()
{
  E0 = parameters["E0"];
  c0 = parameters["c0"];
  W  = parameters["W"];
  r0 = parameters["r0"];
  b  = parameters["b"];
}

// Member variables

double E0;
double c0;
double W;
double r0;
double b;

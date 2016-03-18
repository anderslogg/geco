// Member functions

double ansatz(double E, double L) const
{
  if (E0 <= E)
    return 0.0;

  if (std::abs(L) >= 1.0 / Q)
    return 0.0;
  
  return std::pow(E0 - E, k)*std::pow(1.0 - Q*std::abs(L), l);
}

void init_parameters()
{
  parameters.add("E0", -0.1);
  parameters.add("k",   0.0);
  parameters.add("l",   0.0);
  parameters.add("Q",   0.0);
}

void read_parameters()
{
  E0 = parameters["E0"];
  k  = parameters["k"];
  l  = parameters["l"];
  Q  = parameters["Q"];
}

// Member variables

double E0;
double k;
double l;
double Q;

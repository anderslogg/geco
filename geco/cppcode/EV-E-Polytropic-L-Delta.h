// Member functions

double ansatz(double E, double L) const
{
  if (E >= E0)
    return 0.0;

  if (!_rotation)
    L = std::abs(L);

  if (L <= L0 - l)
    return 0.0; 

  if (L >= L0 + l)
    return 0.0; 

  return std::pow(E0 - E, k);
}

void init_parameters()
{
  parameters.add("E0", 0.9);
  parameters.add("L0", 0.0);
  parameters.add("k",  0.0);
  parameters.add("l",  0.1);
}

void read_parameters()
{
  E0 = parameters["E0"];
  L0 = parameters["L0"];
  k  = parameters["k"];
  l  = parameters["l"];
}

// Member variables

double E0;
double L0;
double k;
double l;

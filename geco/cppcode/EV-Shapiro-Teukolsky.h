// Member functions

double ansatz(double E, double L) const
{
  if (E0 <= E)
    return 0.0;

  if (_rotation && L <= 0.0)
    return 0.0;

  return std::pow( std::pow((E / E0), 2.0), -d)*std::pow( (1 - std::pow(E / E0, 2.0)), d)*(1.0 / (L0*std::sqrt(pi)))*std::exp(- std::pow(L / L0, 2.0));
}

void init_parameters()
{
  parameters.add("E0", 0.9);
  parameters.add("d",  0.0);
  parameters.add("L0", 2.0);
}

void read_parameters()
{
  E0  = parameters["E0"];
  d   = parameters["d"];
  L0  = parameters["L0"];
}

// Member variables

double E0;
double d;
double L0;

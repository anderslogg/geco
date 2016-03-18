// Member functions

double ansatz(double E, double L) const
{
  if (E0 <= E)
    return 0.0;

  if (_rotation && L <= 0.0)
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

// Member variables

double E0;
double k;
double L0;
double pm;

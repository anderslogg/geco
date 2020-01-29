// Member functions

double ansatz(double E, double L) const
{
  if (E0 <= E)
    return 0.0;

  double I = std::pow(E0 - E, k);
  double J = 1.0 - Q*std::abs(L);

  if (J <= 0.0)
    return 0.0;

  else
    return I *= std::pow(J, l);

}

void init_parameters()
{
  parameters.add("E0", -0.1);
  parameters.add("k",   0.0);
  parameters.add("l",   0.0);
  parameters.add("Q",   0.0);
  parameters.add("weight", 1.0);
  parameters.add("model", "VP-E-Polytropic-L-Andreasson");
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
double weight;
std::string model;

// Polytropic ansatz for Vlasov-Poisson

// Ansatz class name
#define CLASSNAME VPEPLP

// Member functions

double ansatz(double E, double L) const
{
  if (E0 <= E)
    return 0.0;

  if (std::abs(L) < L0)
    return 0.0;

  return std::pow(E0 - E, k)*std::pow(std::abs(L) - L0, l);
}

void init_parameters()
{
  parameters.add("E0", -0.1);
  parameters.add("L0",  0.0);
  parameters.add("k",   0.0);
  parameters.add("l",   0.0);
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

// Ansatz for Vlasov-Poisson with Gaussian distribution in L and polytropic in E.

// Member functions

double ansatz(double E, double L) const
{
  if (E0 <= E)
    return 0.0;

  return std::pow(E0 - E, k)*(1.0 / L0)*std::exp(pm*std::pow(L / L0, 2.0));
}

void init_parameters()
{
  parameters.add("E0", -0.1);
  parameters.add("L0",  1.0);
  parameters.add("k",   0.0);
  parameters.add("pm",  -1.0);
  parameters.add("weight", 1.0);
  parameters.add("model", "VP-E-Polytropic-L-Gaussian");
}

void read_parameters()
{
  E0  = parameters["E0"];
  L0  = parameters["L0"];
  k   = parameters["k"];
  pm  = parameters["pm"];

}

// Member variables

double E0;
double k;
double L0;
double pm;
double weight;
std::string model;

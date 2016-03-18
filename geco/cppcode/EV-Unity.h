// Member functions

double ansatz(double E, double L) const
{ 
  return 1.0;
}

void init_parameters()
{
  parameters.add("E0", 2.0);
  parameters.add("L0", 0.0);
  parameters.add("k",  0.0);
  parameters.add("l",  0.0);
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

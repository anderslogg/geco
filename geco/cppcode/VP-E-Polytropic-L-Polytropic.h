// Polytropic ansatz for Vlasov-Poisson

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include "VPAnsatz.h"

class VPEPolyLPoly : public VPAnsatz
{
public:
  // Member functions

  VPEPolyLPoly() : VPAnsatz() 
  {
      // Set default parameter values
      init_parameters();
  };

  double ansatz(double E, double L) const override
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

private:
  // Member variables
  // Note that the energy cutoff E0 is declared within VPAnsatz.

  double L0;
  double k;
  double l;
}; // end class VPEPolyLPoly


PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<VPAnsatz, std::shared_ptr<VPAnsatz>, dolfin::Expression>
    (m, "VPAnsatzTemplate")
    .def("set_fields", (void (VPAnsatz::*)(std::shared_ptr<const dolfin::Function>)) &VPAnsatz::set_fields)
    .def("set_integration_parameters", (void (VPAnsatz::*)(std::size_t)) &VPAnsatz::set_integration_parameters)
    .def("reset", &VPAnsatz::reset)
    .def("radius_of_support", (double (VPAnsatz::*)()) &VPAnsatz::radius_of_support);

  py::class_<VPEPolyLPoly, VPAnsatz, std::shared_ptr<VPEPolyLPoly>>
    (m, "VPAnsatz")
    .def(py::init<>())
    .def("init_parameters", &VPEPolyLPoly::init_parameters)
    .def("read_parameters", &VPEPolyLPoly::read_parameters);
}

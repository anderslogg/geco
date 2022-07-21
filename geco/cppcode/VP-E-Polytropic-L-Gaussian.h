// Ansatz for Vlasov-Poisson with Gaussian distribution in L and polytropic in E.

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include "VPAnsatz.h"

class VPEPolyLGauss : public VPAnsatz
{
public:
  // Member functions

  VPEPolyLGauss() : VPAnsatz() 
  {
      // Set default parameter values
      init_parameters();
  };

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
  }

  void read_parameters()
  {
    E0  = parameters["E0"];
    L0  = parameters["L0"];
    k   = parameters["k"];
    pm  = parameters["pm"];
  }

private:
  // Member variables
  // Note that the energy cutoff E0 is declared within VPAnsatz.

  double k;
  double L0;
  double pm;

}; // end class VPEPolyLGauss


PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<VPAnsatz, std::shared_ptr<VPAnsatz>, dolfin::Expression>
    (m, "VPAnsatzTemplate")
    .def("set_fields", (void (VPAnsatz::*)(std::shared_ptr<const dolfin::Function>)) &VPAnsatz::set_fields)
    .def("set_integration_parameters", (void (VPAnsatz::*)(std::size_t)) &VPAnsatz::set_integration_parameters)
    .def("reset", &VPAnsatz::reset)
    .def("radius_of_support", (double (VPAnsatz::*)()) &VPAnsatz::radius_of_support);

  py::class_<VPEPolyLGauss, VPAnsatz, std::shared_ptr<VPEPolyLGauss>>
    (m, "VPAnsatz")
    .def(py::init<>())
    .def("init_parameters", &VPEPolyLGauss::init_parameters)
    .def("read_parameters", &VPEPolyLGauss::read_parameters);
}
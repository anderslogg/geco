/* Pybind11 code for binding EVAnsatz */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include "EVAnsatz.h"
#include "EV-E-Polytropic-L-Polytropic.h"

PYBIND11_MODULE(SIGNATURE, m)
{
  // EVAnsatz abstract class
  py::class_<EVAnsatz, std::shared_ptr<EVAnsatz>, dolfin::Expression>
    (m, "EVAnsatzTemplate")
    .def("set_fields", (void (EVAnsatz::*)(std::shared_ptr<const dolfin::Function>, std::shared_ptr<const dolfin::Function>, std::shared_ptr<const dolfin::Function>, std::shared_ptr<const dolfin::Function>)) &EVAnsatz::set_fields)
    .def("set_integration_parameters", (void (EVAnsatz::*)(std::size_t)) &EVAnsatz::set_integration_parameters)
    .def("reset", &EVAnsatz::reset)
    .def("radius_of_support", (double (EVAnsatz::*)()) &EVAnsatz::radius_of_support);

  // EV-E-Polytropic-L-Polytropic
  py::class_<EVEPolyLPoly, EVAnsatz, std::shared_ptr<EVEPolyLPoly>>
    (m, "EVEPolytropicLPolytropic")
    .def(py::init<>())
    .def("init_parameters", &EVEPolyLPoly::init_parameters)
    .def("read_parameters", &EVEPolyLPoly::read_parameters);
}

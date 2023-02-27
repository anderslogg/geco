/* Pybind11 code for binding EVAnsatz */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include "EVAnsatz.h"
#include "EV-E-Polytropic-L-Polytropic.h"
#include "EV-E-Polytropic-L-Gaussian.h"
#include "EV-E-Polytropic-L-Andreasson.h"
#include "EV-E-Polytropic-L-Delta.h"
#include "EV-E-Delta-L-Delta.h"

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
 
  // EV-E-Polytropic-L-Gaussian
  py::class_<EVEPolyLGauss, EVAnsatz, std::shared_ptr<EVEPolyLGauss>>
    (m, "EVEPolytropicLGaussian")
    .def(py::init<>())
    .def("init_parameters", &EVEPolyLGauss::init_parameters)
    .def("read_parameters", &EVEPolyLGauss::read_parameters);

  // EV-E-Polytropic-L-Andreasson
  py::class_<EVEPolyLAndreasson, EVAnsatz, std::shared_ptr<EVEPolyLAndreasson>>
    (m, "EVEPolytropicLAndreasson")
    .def(py::init<>())
    .def("init_parameters", &EVEPolyLAndreasson::init_parameters)
    .def("read_parameters", &EVEPolyLAndreasson::read_parameters);

  // EV-E-Polytropic-L-Delta
  py::class_<EVEPolyLDelta, EVAnsatz, std::shared_ptr<EVEPolyLDelta>>
    (m, "EVEPolytropicLDelta")
    .def(py::init<>())
    .def("init_parameters", &EVEPolyLDelta::init_parameters)
    .def("read_parameters", &EVEPolyLDelta::read_parameters);

  // EV-E-Delta-L-Delta
  py::class_<EVEDeltaLDelta, EVAnsatz, std::shared_ptr<EVEDeltaLDelta>>
    (m, "EVEDeltaLDelta")
    .def(py::init<>())
    .def("init_parameters", &EVEDeltaLDelta::init_parameters)
    .def("read_parameters", &EVEDeltaLDelta::read_parameters);
}
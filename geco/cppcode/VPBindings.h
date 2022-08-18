/* Pybind11 code for binding VPAnsatz */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include "VPAnsatz.h"
#include "VP-E-Polytropic-L-Polytropic.h"
#include "VP-E-Polytropic-L-Gaussian.h"
#include "VP-E-Polytropic-L-Andreasson.h"
#include "VP-Evans-L-Polytropic.h"
#include "VP-Rowley.h"

PYBIND11_MODULE(SIGNATURE, m)
{
  // VPAnsatz abstract class
  py::class_<VPAnsatz, std::shared_ptr<VPAnsatz>, dolfin::Expression>
    (m, "VPAnsatzTemplate")
    .def("set_fields", (void (VPAnsatz::*)(std::shared_ptr<const dolfin::Function>)) &VPAnsatz::set_fields)
    .def("set_integration_parameters", (void (VPAnsatz::*)(std::size_t)) &VPAnsatz::set_integration_parameters)
    .def("reset", &VPAnsatz::reset)
    .def("radius_of_support", (double (VPAnsatz::*)()) &VPAnsatz::radius_of_support);

  // VP-E-Polytropic-L-Polytropic
  py::class_<VPEPolyLPoly, VPAnsatz, std::shared_ptr<VPEPolyLPoly>>
    (m, "VPEPolytropicLPolytropic")
    .def(py::init<>())
    .def("init_parameters", &VPEPolyLPoly::init_parameters)
    .def("read_parameters", &VPEPolyLPoly::read_parameters);

  // VP-E-Polytropic-L-Gaussian
  py::class_<VPEPolyLGauss, VPAnsatz, std::shared_ptr<VPEPolyLGauss>>
    (m, "VPEPolytropicLGaussian")
    .def(py::init<>())
    .def("init_parameters", &VPEPolyLGauss::init_parameters)
    .def("read_parameters", &VPEPolyLGauss::read_parameters);  

  // VP-E-Polytropic-L-Andreasson
  py::class_<VPEPolyLAndreasson, VPAnsatz, std::shared_ptr<VPEPolyLAndreasson>>
    (m, "VPEPolytropicLAndreasson")
    .def(py::init<>())
    .def("read_parameters", &VPEPolyLAndreasson::read_parameters);        

  // VP-Evans-L-Polytropic
  py::class_<VPEvansLPoly, VPAnsatz, std::shared_ptr<VPEvansLPoly>>
    (m, "VPEvansLPolytropic")
    .def(py::init<>())
    .def("read_parameters", &VPEvansLPoly::read_parameters);        

  // VP-Rowley
  py::class_<VPRowley, VPAnsatz, std::shared_ptr<VPRowley>>
    (m, "VPRowley")
    .def(py::init<>())
    .def("read_parameters", &VPRowley::read_parameters);        
}

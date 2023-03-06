/* Pybind11 code for binding VPAnsatz */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include "Density2D.h"
#include "Density3D.h"
#include "PointCloud.h"

PYBIND11_MODULE(SIGNATURE, m)
{
  // Density2D 
  py::class_<Density2D, std::shared_ptr<Density2D>, dolfin::Expression>
    (m, "Density2D")
    .def(py::init<>())
    .def("set_density", (void (Density2D::*)(std::shared_ptr<const dolfin::Function>)) &Density2D::set_density);

  // Density3D 
  py::class_<Density3D, std::shared_ptr<Density3D>, dolfin::Expression>
    (m, "Density3D")
    .def(py::init<>())
    .def("set_density", (void (Density3D::*)(std::shared_ptr<const dolfin::Function>)) &Density3D::set_density);
  
  // PointCloud 
  py::class_<PointCloud, std::shared_ptr<PointCloud>, dolfin::Expression>
    (m, "PointCloud")
    .def(py::init<>())
    .def("set_parameters", (void (PointCloud::*)(std::shared_ptr<const dolfin::Function>, double, double, std::size_t, std::size_t)) &PointCloud::set_parameters)
    .def("save_data", (void (PointCloud::*)(std::string)) &PointCloud::save_data);
};
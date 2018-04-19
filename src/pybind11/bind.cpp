#include <pybind11/pybind11.h>
namespace py = pybind11;

void object3D_bind(py::module &);

PYBIND11_MODULE(canspy, m) {
  object3D_bind(m);
}

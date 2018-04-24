#include "cans/object3D.h"

// #include <pcl/point_types.h>
// #include <pcl/point_cloud.h>
// #include <pcl/common/common.h>
// #include <pcl/registration/correspondence_estimation.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

using namespace PLib;

void object3D_bind(py::module &m) {

  py::class_<Object3D>(m, "Object3D")
    .def(py::init<>())
    .def(py::init<int, int,
        Eigen::Array<float,1,Eigen::Dynamic>&,
        Eigen::Array<float,1,Eigen::Dynamic>&,
        Eigen::Array<float,3,Eigen::Dynamic>&,
        int, int>())

    .def("updateObject3D", &Object3D::updateObject3D)
    .def("getDistanceFromPointToSurface", 
        &Object3D::getDistanceFromPointToSurface) 
    // .def("getDistanceFromPointToSurface",[](Object3D& obj, Eigen::Vector3f& query, int ms, int mt){
    //     float dist;
    //     cout << "Distance is " << dist << endl;
    //     cout << "ms " << ms << " mt " << mt << endl;
    //     cout << "query is: " << query << endl;
    //     dist = obj.getDistanceFromPointToSurface(query, ms, mt);
    //     cout << "Distance is " << dist << endl;
    //     return dist;
    // })
    .def("getBatchDistanceFromPointsToSurface", 
        &Object3D::getBatchDistanceFromPointsToSurface);
}
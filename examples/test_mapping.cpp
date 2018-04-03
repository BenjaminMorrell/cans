#include <ros/ros.h>
#include <iostream>

#include <nurbsS.h>
#include <nurbs.h>

#include "cans/mapping3D.h"
#include "cans/object3D.h"

using namespace std;
using namespace PLib;

void test_mapping_class(){

  // Create a matrix "data"
  // Define control points
  Matrix_Point3Df scan(3,3);
  float k;
  
  for (int i = 0; i<scan.rows(); i++){
    for (int j=0; j<scan.cols(); j++){
      k = i*j + 2.0*i - j*3.0;
      scan(i,j) = Point3Df(i,j,k);
    }
  }

  cout << "Data is: " << scan << endl;

  // Initialise class
  Mapping3D map;

  // Get centre
  // Point3Df centre = map.compute_centre_of_data(scan);

  // cout << "Centre is: " << centre << endl;

  // Initialise Object3D
  Object3D obj(scan);
  
  Point3Df centre2 = obj.getCentre();

  cout << "Centre is: " << centre2 << endl;

  Matrix_Point3Df scan2(25,25);
  k = 0;
  
  for (int i = 0; i<scan2.rows(); i++){
    for (int j=0; j<scan2.cols(); j++){
      k = i*j + 2.0*i - j*3.0;
      scan2(i,j) = Point3Df(i,j,k);
    }
  }

  PlNurbsSurfacef surf ;

  // Least Squares Fit
  surf.leastSquares(scan2,3,3,15,15) ;

  // Init object with data
  Object3D obj2(scan2, 3, 3, 15, 15);

  int degu = obj2.degreeU();

  cout << "degree: " << degu << endl;

  float knotu = obj2.knotU(3);

  cout << "knot 3: " << knotu << endl;

  cout << "Point (0.3,0.3) on surface: " << surf.pointAt(0.3,0.3) << "\nPoint (0.3,0.3) on obj: " << obj2.pointAt(0.3,0.3) << endl;

  cout << "Centre from data: " << obj2.getCentre() << endl;

  obj2.computeCentreFromControlPoints();

  cout << "Centre from Control Points: " << obj2.getCentre() << endl;

  cout << "Object Size from data: " << obj2.getObjSize() << endl;
  obj2.computeSizeFromControlPoints();
  cout << "Object Size from Control Points: " << obj2.getObjSize() << endl;

  Object3D obj3; // don't empty declare with open and closed brackets!

  cout << "Default centre: " << obj3.getCentre() << endl;

  obj3.leastSquares(scan2,3,3,15,15) ;

  cout << "Point (0.3,0.3) on obj3: " << obj3.pointAt(0.3,0.3) << endl;

}

// void testMappingWorkflow(int argc, char** argv){

//   pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
//   pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
//   pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>);
//   pcl::PointCloud<pcl::PointNormal>::Ptr cloud2 (new pcl::PointCloud<pcl::PointNormal>);
 
 
//   // Fill in the cloud data
//   pcl::PCDReader reader;
//   reader.read (argv[1], *cloud_blob);
//   if (argc == 3){
//     reader.read (argv[2], *cloud_blob2);
//   }
  

//   // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
//   pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
//   if (argc > 3){
//     pcl::fromPCLPointCloud2 (*cloud_blob, *cloud2_orig); 
//   }else{
//     pcl::fromPCLPointCloud2 (*cloud_blob2, *cloud2); 
//   }

// }


int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "test_mapping");
  ros::NodeHandle nh;

  // FIt a NURBS surface
  test_mapping_class();

  // Spin
//   ros::spin ();
}

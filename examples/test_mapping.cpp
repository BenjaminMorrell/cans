#include <ros/ros.h>
#include <iostream>

#include <nurbsS.h>
#include <nurbs.h>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include <pcl/io/pcd_io.h>

#include "cans/mapping3D.h"
#include "cans/object3D.h"
#include "cans/splitSurface.h"

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

void testMappingWorkflow(int argc, char** argv){

  if (argc < 3){
    cout << "Incorrect input: need to input filepaths to two point clouds" << endl;
    return;
  }

  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>); 
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2 (new pcl::PointCloud<pcl::PointXYZ>);

 
 
  // Fill in the cloud data
  pcl::PCDReader reader;
  reader.read (argv[1], *cloud_blob);
  reader.read (argv[2], *cloud_blob2);
  
  // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
  pcl::fromPCLPointCloud2 (*cloud_blob2, *cloud2); 

  // Init Mapping class
  Mapping3D mp;
  mp.numRowsDesired = 45;
  mp.numColsDesired = 45;
  mp.maxNanAllowed = 10;

  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numRowsDesired));
  
  // Process Scan 1
  mp.meshFromScan(cloudReduced, cloud);

  // Get Mesh in format for NURBS
  Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(cloudReduced);

  // Add object
  Object3D obj(mesh);

  cout << "Point (0.3,0.3) on obj: " << obj.pointAt(0.3,0.3) << endl;
  
  // Write NURBS
  obj.writeVRML("meshTest.wrl",Color(255,100,255),50,80);

  // Mesh from scan2
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced2 (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numRowsDesired));
  mp.meshFromScan(cloudReduced2, cloud2);

  // Split new surface observation
  SplitSurface ss;

  ss.splitNewSurfaceObservation(cloudReduced, cloudReduced2);

  cout << "Extend Direction is: " << ss.extendDirection << endl;

  // Get matrix from new data
  Matrix_Point3Df mesh2 = mp.nurbsDataFromPointCloud(cloudReduced2, ss.newDataIndices);

  cout << "Mesh2 size is (" << mesh2.rows() << ", " << mesh2.cols() << ")\n";

  cout << "Mesh2 is: " << mesh2 << endl;

  std::vector<int> nCtrlNew = mp.computeNumberOfControlPoints(ss.extendDirection, mesh2, obj);

  cout << "new control points: " << nCtrlNew[0] << ", " << nCtrlNew[1] << endl;

  // Create Object for new surface? (nknots = deg + nctrl + 1)
  Object3D obj2(mesh2, 3, 3, nCtrlNew[0], nCtrlNew[1]);
  obj2.writeVRML("mesh2Test.wrl",Color(255,100,255),50,80);

  cout << "Point (0.3,0.3) on obj2: " << obj2.pointAt(0.3,0.3) << endl;
  cout << "Point (0.6,0.) on obj2: " << obj2.pointAt(0.6,0.0) << endl;
  cout << "Point (0.,0.) on obj2: " << obj2.pointAt(0.0,0.0) << endl;

  // Join Surfaces
  Object3D obj3 = mp.joinSurfaces(obj,obj2,ss.extendDirection);

  cout << "Point (0.3,0.3) on obj3: " << obj3.pointAt(0.3,0.3) << endl;

  // obj3.writeVRML("mesh3Test.wrl",Color(255,100,255),50,80);

  mp.updateObject(obj, cloudReduced, cloudReduced2);

}

void testDataAssociation(int argc, char ** argv){
  if (argc < 3){
    cout << "Incorrect input: need to input filepaths to two point clouds" << endl;
    return;
  }

  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>); 
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2 (new pcl::PointCloud<pcl::PointXYZ>);

 
 
  // Fill in the cloud data
  pcl::PCDReader reader;
  reader.read (argv[1], *cloud_blob);
  reader.read (argv[2], *cloud_blob2);
  
  // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
  pcl::fromPCLPointCloud2 (*cloud_blob2, *cloud2); 

  // Create objects
  // Init Mapping class
  Mapping3D mp;
  mp.numRowsDesired = 45;
  mp.numColsDesired = 45;
  mp.maxNanAllowed = 10;

  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numRowsDesired));
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced2 (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numRowsDesired));
  
  // Process Scan 1
  mp.meshFromScan(cloudReduced, cloud);
  mp.meshFromScan(cloudReduced2, cloud2);

  // // Get Mesh in format for NURBS
  // Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(cloudReduced);
  // Matrix_Point3Df mesh2 = mp.nurbsDataFromPointCloud(cloudReduced2);

  // // Create Objects and add to map
  // mp.objectMap.push_back(Object3D(mesh));
  // mp.objectMap.push_back(Object3D(mesh2));

  float arr[] = {0.1,0.2,0.4,5.0,0,0,0};
  vector<float> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );

  mp.addObject(cloudReduced,vec);

  vec[0] = 1.2;
  vec[1] = 0.2;
  vec[2] = 2.2;
  vec[3] = 1.2;

  mp.addObject(cloudReduced2, vec);

  cout << "Object 1 size: " << mp.objectMetrics[0][3] << endl;
  cout << "Object 1 size: " << mp.objectMap[0].getObjSize() << endl;
  cout << "Object 2 size: " << mp.objectMetrics[1][3] << endl;

  // cout << "Object 1 at (0,0): " << mp.objectMap[0]->pointAt(0.0,0.0) << endl;
  // cout << "Object 2 at (0,0): " << mp.objectMap[1].pointAt(0.0,0.0) << endl;

  Matrix_HPoint3Df controlPoints = mp.objectMap[0].ctrlPnts();

  cout << "Control point (0,0): " << controlPoints(0,0) << endl;
  cout << "Control point (1,1): " << controlPoints(1,1) << endl;
  cout << "Control point rows: " << controlPoints.rows() << endl;
  cout << "Control point cols: " << controlPoints.rows() << endl;

  cout << "Obj0 order: " << mp.objectMap[0].degreeU() << endl;
  cout << "Obj0 knots: " << mp.objectMap[0].knotU() << endl;

  int id = mp.dataAssociation(mp.objectMetrics[1]);

  mp.dataAssociation(mp.objectMetrics[0]);


  // std::vector<std::vector<float> > bla;

  // // vector<float> vec1 = {1,2,3,4,5,6,7};

  // bla.push_back(vec);

  // bla.push_back(vec);

  
  // cout << "Object metrics 0, 0: " << bla[0][0] << endl;
  // cout << "Object metrics 1, 5: " << bla[1][5] << endl;
  // cout << "Size dimension 1 is: " << bla.size() << endl;
  // cout << "Size dimension 2 is: " << bla[0].size() << endl;

  // std::vector<Object3D> objectMap;

  // objectMap.push_back(Object3D());
  // objectMap.push_back(Object3D());

  // cout << "Object map 0: " << objectMap[0].getCentre() << endl;
  // cout << "Object map 1: " << objectMap[1].getCentre() << endl;
  // cout << "Number of objects is: " << objectMap.size() << endl;
}


int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "test_mapping");
  ros::NodeHandle nh;

  // FIt a NURBS surface
  // test_mapping_class();

  // testMappingWorkflow(argc, argv);

  testDataAssociation(argc,argv);

  // Spin
//   ros::spin ();
}

#include <ros/ros.h>
#include <iostream>
#include <fstream>
#include <iomanip>

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
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2 (new pcl::PointCloud<pcl::PointNormal>);

 
 
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

  // mp.updateObject(obj, cloudReduced, cloudReduced2);

}

void testDataAssociation(int argc, char ** argv){
  if (argc < 3){
    cout << "Incorrect input: need to input filepaths to two point clouds" << endl;
    return;
  }

  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2 (new pcl::PointCloud<pcl::PointNormal>);

 
 
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

  std::vector<float> searchMetrics = mp.computeSearchMetrics(cloudReduced);

  cout << "search metrics are:\n";
  for (int i = 0; i < 7; i ++){cout << searchMetrics[i] << "\t";}cout << endl;

  id = mp.dataAssociation(searchMetrics);

  cout << "Cloud compute metrics search gives: " << id << endl;

  id = mp.dataAssociation(searchMetrics);

  cout << "Cloud compute metrics search gives (should be 0): " << id << endl;

  // Test with object
  // Get Mesh in format for NURBS
  Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(cloudReduced2);

  // Add object
  Object3D obj(mesh);

  // Compute metrics
  searchMetrics = mp.computeSearchMetrics(obj);

  id = mp.dataAssociation(searchMetrics);

  cout << "search metrics are:\n";
  for (int i = 0; i < 7; i ++){cout << searchMetrics[i] << "\t";}cout << endl;

  cout << "Cloud compute metrics search gives (should be 1): " << id << endl;


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




// void processInputPCs(int argc, char ** argv, pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloud2){
  
  
//   cout << "Loaded clouds" << endl;

// }

void testHighLevelFunctions(int argc, char ** argv){

  if (argc < 3){
    cout << "Incorrect input: need to input filepaths to two point clouds" << endl;
    return;
  }
  
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2 (new pcl::PointCloud<pcl::PointNormal>);
   
  // Fill in the cloud data
  pcl::PCDReader reader;
  reader.read (argv[1], *cloud_blob);
  reader.read (argv[2], *cloud_blob2);
  
  // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
  pcl::fromPCLPointCloud2 (*cloud_blob2, *cloud2); 

  

  // processInputPCs(argc, argv, cloud, cloud2);

    // Create objects
  // Init Mapping class
  Mapping3D mp;
  mp.numRowsDesired = 45;
  mp.numColsDesired = 45;
  mp.maxNanAllowed = 10;
  
  // Transform
  Eigen::Affine3f transform = Eigen::Affine3f::Identity();
  // transform.translation() << atof(argv[2]), atof(argv[3]), atof(argv[4]);
    
  // // Rotation
  // // Rot vec
  // Eigen::Vector3f rotVec(atof(argv[6]),atof(argv[7]),atof(argv[8]));
  // rotVec.normalize();

  // transform.rotate (Eigen::AngleAxisf (atof(argv[5]), rotVec));

  int res = mp.processScan(cloud,transform);

  // Test if the object exists
  cout << "Point (0,0) on object 1 is: " << mp.objectMap[0].pointAt(0.0,0.0) << endl;
  cout << "Point (0.4,0.4) on object 1 is: " << mp.objectMap[0].pointAt(0.4,0.4) << endl;

  cout << "Object 0 knotU: " << mp.objectMap[0].knotU() << endl;

  // Second scan
  res = mp.processScan(cloud2, transform);

  cout << "res from process scan is: " << res << endl;


  
}

Eigen::Array<float,6,10> readPathTextFile(char * filename){

  std::ifstream file;
  std::string line;

  Eigen::Array<float,6,10> state(6,10);

  float f1, f2, f3, f4, f5, f6, f7;
  char* s1[10], s2[7], s3[6], s4[2], s5[2], s6[2], s7[3], s8[3], s9[3];

  FILE * pFile;

  pFile = fopen (filename, "r");

  for (int i = 0; i < 10; i++){
    fscanf(pFile, "%s%f%s%s%s%f%s%f%s%f%s%f%s%f%s%f", &s1, &f1, &s2, &s3, &s4, &f2, &s5, &f3, &s6, &f4, &s7, &f5, &s8, &f6, &s9, &f7);

    // cout << f1 << ", " << f2 << ", " << f3 << ", " << f4 << ", " << f5 << endl;
    // cout << s1 << endl;
    // cout << s2 << endl;
    // cout << s3 << endl;
    // cout << s4 << endl;
    // cout << s5 << endl;
    // cout << s6 << endl;
    // cout << s7 << endl;
    // cout << s8 << endl;
    // cout << s9 << endl;
    state(0,i) = f2;
    state(1,i) = f3;
    state(2,i) = f4;
    state(3,i) = f5;
    state(4,i) = f6;
    state(5,i) = f7;
    
  }
  

  // cout << f1 << ", " << f2 << ", " << f3 << ", " << f4 << ", " << f5 << endl;

  fclose(pFile);

  cout << state << endl;

  return state;


}

void testBlenderSequence(int argc, char ** argv){
  
  int numberOfScans = 1;
  std::string filename;
  std::string filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Data000";
  std::string outFilename;
  std::string outFilestem = "blob_scan_res_new_";
  int res;
  Eigen::Affine3f transform = Eigen::Affine3f::Identity();
  // Eigen::Matrix3f rotMat;

  cout << "In blender sequence test " << endl;

  // Load state
  Eigen::Array<float,6,10> state = readPathTextFile("/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Path.txt");
  cout << "State is:\n" << state << endl;

  // Init Mapping class
  Mapping3D mp;
  mp.numRowsDesired = 95;
  mp.numColsDesired = 95;
  mp.maxNanAllowed = 10;
  mp.removeNanBuffer = 3;
  mp.nCtrlDefault[0] = 30;
  mp.nCtrlDefault[1] = 30;

  mp.newRowColBuffer = 20; // How many non new points in a row or column are permissible

  pcl::PCDReader reader;

  if (argc > 1){  
    numberOfScans = atoi(argv[1]);
  }

  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 

  for (int i = 0; i < numberOfScans; i++){
    cout << "Processing Scan " << i << endl;
    cout << "State x is: " << state(0,i) << endl;

    // Get filename:
    if (i < 9){
      filename = filestem + "0" + static_cast<ostringstream*>( &(ostringstream() << (i+1)) )->str() + ".pcd";
      cout << "filename is : " << filename;
    }else{
      filename = filestem + static_cast<ostringstream*>( &(ostringstream() << (i+1)) )->str() + ".pcd";
      cout << "filename is : " << filename;
    }

    // Read Scan
    reader.read (filename, *cloud_blob);
    cout << "scan read" << endl;

    // Convert to PCL cloud
    pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 

    // cout << "Cloud blob transform is: " << cloud_blob->sensor_origin_ << "\n" << cloud_blob->sensor_orientation_ << endl;

    // Fill transform
    transform = Eigen::Affine3f::Identity();

    // Translation
    transform(0,3) = state(0,i);
    transform(1,3) = state(1,i);
    transform(2,3) = state(2,i);
    
    cout << "transform is " << transform.matrix() << endl;

    // 3, 2, 1 Euler transformation
    transform.rotate (Eigen::AngleAxisf(state(5,i),Eigen::Vector3f::UnitZ()));
    transform.rotate (Eigen::AngleAxisf(state(4,i),Eigen::Vector3f::UnitY()));
    transform.rotate (Eigen::AngleAxisf(state(3,i),Eigen::Vector3f::UnitX()));

    cout << "transform is " << transform.matrix() << endl;

    // Process Scan
    res = mp.processScan(cloud,transform);

    // Save Scan
    outFilename = outFilestem + static_cast<ostringstream*>( &(ostringstream() << (i+1)) )->str() + ".wrl";
    mp.objectMap[0].writeVRML(outFilename.c_str(),Color(255,100,255),50,80);  

    


  }

  cout << "\nSize of Object map is: " << mp.objectMap.size() << endl;
  
  // Write result to pcd
  mp.writeObjectPCDFile("endResultBlob.pcd", 0, 125, 125);


  if (argc > 2){
    mp.objectMap[0].write(argv[2]);    
    cout << "Starting centre is: " << mp.objectMap[0].getCentre() << ", and centre: " << mp.objectMap[0].getObjSize() << endl;

    Object3D obj;
    obj.readObject3D(argv[2]);
    // obj.read("testSaveObject");
    // obj.computeSizeFromControlPoints();
    // obj.computeCentreFromControlPoints(); // TODO - color

    obj.writeVRML("savedObject.wrl",Color(255,100,255),50,80); 
    cout << "Loaded object has centre: " << obj.getCentre() << ", color: " << obj.getColor() << ", size: " << obj.getObjSize() << endl;
  }
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

  // testDataAssociation(argc,argv);

  // testHighLevelFunctions(argc, argv);

  testBlenderSequence(argc, argv);

  // readPathTextFile("/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Path.txt");

  // Spin
//   ros::spin ();
}

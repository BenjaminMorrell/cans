#include <ros/ros.h>
#include <iostream>

#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <vector>
#include <ctime>
#include <cmath>

#include "cans/splitSurface.h"

using namespace std;


void splitSurfaceRun(int argc, char ** argv){

  if (argc < 3){
    cout << "Incorrect Usage: Need to input two files for map and observation clouds" << endl;
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

  // test loading in object
  SplitSurface ss;

  ss.setInputMap(cloud);
  ss.setInputObservation(cloud2);
  cout << "Cloud width in split surface is: " << ss.getCloudWidth() << endl;

  ss.findNonOverlappingData();
  cout << "Number of new points: " << ss.newPointsArray.count() << endl;

  // char * something = "hello";

  std::string stng = "jello";

  stng[1] = 'L';

  cout << stng << endl;

  ss.getNewDataExtendDirection();

  cout << "Extend Direction is: " << ss.extendDirection << endl;

  cout << "Indices for new data are: " << ss.newDataIndices << endl;
  cout << "Indices for overlap data are: " << ss.overlapDataIndices << endl;

  Eigen::Vector3f vec;
  Eigen::Vector3f vec2;

  vec = cloud->at(0,0).getArray3fMap();
  vec2 = cloud2->at(0,0).getArray3fMap() - cloud2->at(0,1).getArray3fMap();

  cout << "vector Eigen is: " << vec << endl;
  cout << "vector2 Eigen is: " << vec2 << endl;

  cout << "Dot product is: " << vec.dot(vec2) << endl;

  cout << "Norms are: 1: " << vec.norm() << ", 2: " << vec2.norm() << endl;

  float theta = acos(vec.dot(vec2)/vec.norm()/vec2.norm());

  cout << "Theta is: " << theta << endl;

  ss.getMapDataExtendDirection();

  cout << "Extend Direction is: " << ss.extendDirection << endl;

  SplitSurface ss2;

  ss2.splitNewSurfaceObservation(cloud, cloud2);

  cout << "Extend Direction is: " << ss2.extendDirection << endl;

  // Create data
  // pcl::PointCloud<pcl::PointNormal>::Ptr pcTest (new pcl::PointCloud<pcl::PointNormal>(10,10));
  // int n_points = 10;
  // for (int i = 0; i<n_points; i++){
  //     for (int j = 0; j<n_points; j++){
  //         pcTest->at(j,i).x = (float)j/(n_points-1);
  //         pcTest->at(j,i).y = (float)i/(n_points-1);
  //     } 
  // }
  // cout << "pcTest at (5,4): " << pcTest->at(5,4) << endl;

  // Test new function to 

}



int main (int argc, char ** argv){

  ros::init (argc, argv, "splitSurfaceTest");
  ros::NodeHandle nh;

  splitSurfaceRun(argc,argv);
  // kdtreeSearch();
}
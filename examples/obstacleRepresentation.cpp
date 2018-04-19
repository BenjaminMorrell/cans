#include <ros/ros.h>
#include <iostream>

//PCL includes
#include <Eigen/Core>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/time.h>
#include <pcl/io/pcd_io.h>

#include <pcl/common/transforms.h>


#include <cans_msgs/Object3D.h>

#include "cans/mapping3D.h"
#include "cans/object3D.h"
#include "cans/splitSurface.h"

using namespace std;


void testGetDistancesToSurface(int argc, char ** argv){

  // Init mapping class
  Mapping3D mp;
  mp.numRowsDesired = 45;
  mp.numColsDesired = 45;
  mp.maxNanAllowed = 10;
  mp.removeNanBuffer = 2;

  int msSurf = 55;
  int mtSurf = 155;

  // Load the object
  mp.addObjectFromFile("blob_object.obj");

  // Get a point cloud
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>(mtSurf,msSurf,pcl::PointNormal())); 

  mp.objectMap[0].getSurfacePointCloud(cloud,msSurf,mtSurf);

  // Write to pcd
  pcl::PCDWriter writer;
  writer.write<pcl::PointNormal> ("pcl_from_Object3D.pcd", *cloud, false);

  // Get distance to a point 
  Eigen::Vector3f query; 
  query << 0.5,-1.0,0.65;
  cout << "query is " << query << endl;

  float signedDist = mp.objectMap[0].getDistanceFromPointToSurface(query);

  cout << "Signed distance out is: " << signedDist << endl;

  //----------------------------
  // Test batch of points
  int nPoints = 35;
  Eigen::Array<float,3,Eigen::Dynamic> queryBatch(3,nPoints);
  Eigen::Array<float,1,Eigen::Dynamic> distBatch(1,nPoints);

  float ratio;
  for (int i = 0; i < nPoints; i++){
    ratio = (float)i/(float)nPoints;
    queryBatch(0,i) = -0.5 + 1.8*ratio;
    queryBatch(1,i) = -1.2 + 1.6*ratio;
    queryBatch(2,i) = 0.5;
  }

  distBatch = mp.objectMap[0].getBatchDistanceFromPointsToSurface(queryBatch,101,101);

  cout << "Batch of distances is: " << distBatch << endl;

  for (int i = 0; i < nPoints; i++){
    ratio = (float)i/(float)nPoints;
    queryBatch(0,i) = 0.4;
    queryBatch(1,i) = -0.6;
    queryBatch(2,i) = -0.1 + 1.3*ratio;
  }

  distBatch = mp.objectMap[0].getBatchDistanceFromPointsToSurface(queryBatch,101,101);

  cout << "Batch of distances 2 is: " << distBatch << endl;
}


cans_msgs::Object3D fillObject3DMessage(){

  cans_msgs::Object3D msg;

  msg.ID = 1;

  msg.degU = 3;
  msg.degV = 3;
  msg.nCtrlS = 7;
  msg.nCtrlT = 7;

  Eigen::Array<float,1,11> knots(1,11);
  knots << 0,0,0,0,0.1,0.3,0.7,1,1,1,1;

  // Set size of knotU
  static const float arr[] = {0,0,0,0,0.1,0.3,0.7,1,1,1,1};
  std::vector<float> data (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  
  msg.knotU = data;
  msg.knotV = data;

  // Eigen::Array<float,1,Eigen::Dynamic> controlX(1,49);
  // Eigen::Array<float,1,Eigen::Dynamic> controlY(1,49);
  // Eigen::Array<float,1,Eigen::Dynamic> controlZ(1,49);

  std::vector<float> controlX;
  std::vector<float> controlY;
  std::vector<float> controlZ;

  int index = 0;
  for (int i = 0; i < 7; i++){
    for (int j = 0; j < 7; j++){
      controlX.push_back((float)j);
      controlY.push_back((float)i);
      controlZ.push_back(0.0);
    }
  }

  msg.ControlPoints_x = controlX;
  msg.ControlPoints_y = controlY;
  msg.ControlPoints_z = controlZ;
  
  return msg;
}


int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "planar_seg_testing");
  ros::NodeHandle nh;

  testGetDistancesToSurface(argc,argv);

  ros::Publisher pub = nh.advertise<cans_msgs::Object3D>("object",1,false);

  cans_msgs::Object3D msg = fillObject3DMessage();

  ros::Rate r(0.5);
  while (nh.ok()){
    pub.publish(msg);
    ros::spinOnce();
    r.sleep();
  }
  

}
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
  Eigen::Array<float,4,Eigen::Dynamic> distBatch(4,nPoints);

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

  // Eigen::Array<float,1,11> knots(1,11);
  // knots << 0,0,0,0,0.1,0.3,0.7,1,1,1,1;

  // Set size of knotU
  static const float arr[] = {0,0,0,0,0.1,0.3,0.7,1,1,1,1};
  // static const float arr[] = {0,0,0,0,0.1,0.2,0.3,0.4,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.7,0.76,0.8,0.86,0.9,0.97,1,1,1,1};
  // std::vector<float> data (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  std::vector<float> data;
  for (int i = 0; i < 11 ; i++){
    data.push_back(arr[i]);
  }
  
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

cans_msgs::Object3D fillObject3DMessage2(){

  cans_msgs::Object3D msg;

  msg.ID = 1;

  msg.degU = 3;
  msg.degV = 3;
  msg.nCtrlS = 30;
  msg.nCtrlT = 75;

  // Eigen::Array<float,1,11> knots(1,11);
  // knots << 0,0,0,0,0.1,0.3,0.7,1,1,1,1;

  // Set size of knotU
  // static const float arr[] = {0,0,0,0,0.1,0.3,0.7,1,1,1,1};
  // static const float arr[] = {0,0,0,0,0.1,0.2,0.3,0.4,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.7,0.76,0.8,0.86,0.9,0.97,1,1,1,1};
  // std::vector<float> data (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  std::vector<float> knotU;
  float step = 1.0/(msg.nCtrlS - 1 - msg.degU + 1);
  float val = step;
  for (int i = 0; i < msg.nCtrlS + 1 + msg.degU ; i++){
    if (i < 4){
      knotU.push_back(0.0);
    }else if (i >= msg.nCtrlS){
      knotU.push_back(1.0);
    }else{
      knotU.push_back(val);
      val = val + step;
    }
  }

  std::vector<float> knotV;
  step = 1.0/(msg.nCtrlT - 1 - msg.degV + 1);
  val = step;
  for (int i = 0; i < msg.nCtrlT + 1 + msg.degV ; i++){
    if (i < 4){
      knotV.push_back(0.0);
    }else if (i >= msg.nCtrlT){
      knotV.push_back(1.0);
    }else{
      knotV.push_back(val);
      val = val + step;
    }
    cout << "KnotV is at " << i << " is " << knotV[i] << endl;
  }
  
  msg.knotU = knotU;
  msg.knotV = knotV;

  

  // Eigen::Array<float,1,Eigen::Dynamic> controlX(1,49);
  // Eigen::Array<float,1,Eigen::Dynamic> controlY(1,49);
  // Eigen::Array<float,1,Eigen::Dynamic> controlZ(1,49);

  std::vector<float> controlX;
  std::vector<float> controlY;
  std::vector<float> controlZ;

  int index = 0;
  for (int i = 0; i < msg.nCtrlS; i++){
    for (int j = 0; j < msg.nCtrlT; j++){
      controlX.push_back((float)j/(float)msg.nCtrlT);
      controlY.push_back((float)i/(float)msg.nCtrlS);
      controlZ.push_back(0.0);
    }
  }

  msg.ControlPoints_x = controlX;
  msg.ControlPoints_y = controlY;
  msg.ControlPoints_z = controlZ;


  // Test it can create a valid object
  Object3D obj;

  obj.updateObject3DCPP(msg.degU,msg.degV,knotU,knotV,controlX,controlY,controlZ,msg.nCtrlS,msg.nCtrlT);

  cout << "updated object" << endl;

  cout << "Surface point is " << obj.pointAt(0.0,0.0) << endl;
  

  return msg;
}

cans_msgs::Object3D fillObject3DMessageWithBlob(){

  // Init mapping class
  Mapping3D mp;
  mp.numRowsDesired = 45;
  mp.numColsDesired = 45;
  mp.maxNanAllowed = 10;
  mp.removeNanBuffer = 2;

  int msSurf = 55;
  int mtSurf = 155;

  // Load the object
  mp.addObjectFromFile("/home/bjm/SpaceCRAFT/ros_ws/blob_object.obj");

  cans_msgs::Object3D msg;

  msg.ID = 0;

  msg.degU = mp.objectMap[0].degreeU();
  msg.degV = mp.objectMap[0].degreeV();
  msg.nCtrlS = mp.objectMap[0].ctrlPnts().rows();
  msg.nCtrlT = mp.objectMap[0].ctrlPnts().cols();

  std::vector<float> knotU;
  std::vector<float> knotV;
  for (int i = 0; i < mp.objectMap[0].knotU().size(); i++){
    knotU.push_back(mp.objectMap[0].knotU()[i]);
  }
  for (int i = 0; i < mp.objectMap[0].knotV().size(); i++){
    knotV.push_back(mp.objectMap[0].knotV()[i]);
    cout << "KnotV is at " << i << " is " << knotV[i] << endl;
  }

  cout << "Have knots" << endl;
  // cout << "KnotV is:\n" << knotV << endl;
  
  msg.knotU = knotU;
  msg.knotV = knotV;

  std::vector<float> controlX;
  std::vector<float> controlY;
  std::vector<float> controlZ;

  for (int i = 0; i < msg.nCtrlS; i++){
    for (int j = 0; j < msg.nCtrlT; j++){
      controlX.push_back(mp.objectMap[0].ctrlPnts()(i,j).x());
      controlY.push_back(mp.objectMap[0].ctrlPnts()(i,j).y());
      controlZ.push_back(mp.objectMap[0].ctrlPnts()(i,j).z());
      // controlX.push_back((float)i);
      // controlY.push_back((float)j);
      // controlZ.push_back(0.0);
    }
  }

  cout << "have control points" << endl;

  msg.ControlPoints_x = controlX;
  msg.ControlPoints_y = controlY;
  msg.ControlPoints_z = controlZ;

  Object3D obj;

  obj.updateObject3DCPP(msg.degU,msg.degV,knotU,knotV,controlX,controlY,controlZ,msg.nCtrlS,msg.nCtrlT);

  cout << "updated object" << endl;

  cout << "Surface point is " << obj.pointAt(0.0,0.0) << endl;  
  
  return msg;
}

void messageCallback(cans_msgs::Object3D msg){


  int degU = msg.degU;
  int degV = msg.degV;
  int nCtrlS = msg.nCtrlS;
  int nCtrlT = msg.nCtrlT;

  std::vector<float> knotU = msg.knotU;
  std::vector<float> knotV = msg.knotV;

  std::vector<float> controlX = msg.ControlPoints_x;
  std::vector<float> controlY = msg.ControlPoints_y;
  std::vector<float> controlZ = msg.ControlPoints_z;

  cout << "read message" << endl;

  Object3D obj;


  obj.updateObject3DCPP(degU,degV,knotU,knotV,controlX,controlY,controlZ,nCtrlS,nCtrlT);

  cout << "updated object" << endl;

  cout << "Surface point is " << obj.pointAt(0.0,0.0) << endl;

  // Number of points
  int msSurf = 25;
  int mtSurf = 25;

  // Init a point cloud
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>(mtSurf,msSurf,pcl::PointNormal())); 

  // Get the point cloud
  obj.getSurfacePointCloud(cloud,msSurf,mtSurf);

  cout << "Have point cloud, one point is: " << cloud->at(0,0) << endl;
  

}


int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "planar_seg_testing");
  ros::NodeHandle nh;

  // testGetDistancesToSurface(argc,argv);

  ros::Publisher pub = nh.advertise<cans_msgs::Object3D>("object",1,false);

  // ros::Subscriber sub = nh.subscribe("object", 1, messageCallback);

  // cans_msgs::Object3D msg = fillObject3DMessage2();
  cans_msgs::Object3D msg = fillObject3DMessageWithBlob();


  ros::Rate r(0.2);
  while (nh.ok()){
    pub.publish(msg);
    ros::spinOnce();
    r.sleep();
  }
  

}
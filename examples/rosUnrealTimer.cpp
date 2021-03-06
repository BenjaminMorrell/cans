
#include <ros/ros.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <nurbsS.h>
#include <nurbs.h>

// PCL specific includes
#include <sensor_msgs/PointCloud2.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

#include <tf/transform_listener.h>
#include "tf/transform_datatypes.h"
#include "tf_conversions/tf_eigen.h"
#include "eigen_conversions/eigen_msg.h"

// CANS includes
#include "cans/mapping3D.h"
#include "cans/object3D.h"
#include "cans/splitSurface.h"

using namespace std;
using namespace PLib;

// Init Mapping class
Mapping3D mp;

bool runFlag;
int scanNumber;

pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 

Eigen::Affine3d transformAff;

void 
cloud_cb (const sensor_msgs::PointCloud2ConstPtr& cloud_msg)
{

  if (runFlag){
    // Set flag to false until finished
    runFlag = false;
    cout << "\nProcessing Scan " << scanNumber << "./n/n";
  
    tf::TransformListener listener;
    tf::StampedTransform transform;

    try{
      // listener.waitForTransform("body", "world", cloud_msg->header.stamp, ros::Duration(0.1));
      // listener.lookupTransform("body", "world",  
      //                           cloud_msg->header.stamp, transform);
      listener.waitForTransform("Cam_optical", "world", ros::Time(0), ros::Duration(0.5));
      listener.lookupTransform("Cam_optical", "world",  
                                ros::Time(0), transform);                          
    }
    catch (tf::TransformException ex){
      ROS_ERROR("%s",ex.what());
      // ros::Duration(1.0).sleep();
    }
    // Fill transform
    tf::transformTFToEigen(transform,transformAff);

    cout << "transform is " << transformAff.matrix() << endl;
    
    // Container for original & filtered data
    pcl::PCLPointCloud2* cloud_blob = new pcl::PCLPointCloud2; 
    
    
    // Convert to PCL data type
    pcl_conversions::toPCL(*cloud_msg, *cloud_blob);

    // Convert to Point Cloud <T>
    pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
    
    std::cerr << "PointCloud Received: size: " << cloud->width * cloud->height << endl;
    
    // Set flag back to true - ready for next frame
    runFlag = true;
    
  }
}

// Process Scan
void processPointCloud(const ros::TimerEvent&){
  int res;
  res = mp.processScan(cloud, transformAff.cast<float>());

  std::string filename = "testNURBS_Unreal_" + static_cast<ostringstream*>( &(ostringstream() << (scanNumber)) )->str() + ".wrl";
  mp.objectMap[mp.objectMap.size()-1].writeVRML(filename.c_str(),Color(255,100,255),50,80); 

  cout << "\nFinished Processing Scan " << scanNumber << ".\n\n";
  scanNumber ++;
}

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "NURBS_CANS");
  ros::NodeHandle nh;

  mp.numRowsDesired = 45;
  mp.numColsDesired = 45;
  mp.maxNanAllowed = 10;
  mp.removeNanBuffer = 3;
  mp.nCtrlDefault[0] = 20;
  mp.nCtrlDefault[1] = 20;

  runFlag = true;
  scanNumber = 1;

  // Create a ROS subscriber for the input point cloud
  ros::Subscriber sub = nh.subscribe ("/camera/points2", 1, cloud_cb);

  ros::Timer timer = nh.createTimer(ros::Duration(10.0), processPointCloud);
  
  // Spin
  ros::spin ();
}
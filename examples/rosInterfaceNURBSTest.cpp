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

void 
cloud_cb (const sensor_msgs::PointCloud2ConstPtr& cloud_msg)
{
  
  tf::TransformListener listener;
  tf::StampedTransform transform;

  try{
    listener.waitForTransform("body", "world", cloud_msg->header.stamp, ros::Duration(0.5));
    listener.lookupTransform("body", "world",  
                              cloud_msg->header.stamp, transform);
  }
  catch (tf::TransformException ex){
    ROS_ERROR("%s",ex.what());
    // ros::Duration(1.0).sleep();
  }
  // Fill transform
  Eigen::Affine3d transformAff;

  tf::transformTFToEigen(transform,transformAff);

  cout << "transform is " << transformAff.matrix() << endl;
  

  // Container for original & filtered data
  pcl::PCLPointCloud2* cloud_blob = new pcl::PCLPointCloud2; 
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  
  // Convert to PCL data type
  pcl_conversions::toPCL(*cloud_msg, *cloud_blob);

  // Convert to Point Cloud <T>
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
   
  std::cerr << "PointCloud Received: size: " << cloud->width * cloud->height << endl;

  

  // Process Scan
  int res;
  res = mp.processScan(cloud, transformAff.cast<float>());

  mp.objectMap[0].writeVRML("testNURBSROS.wrl",Color(255,100,255),50,80); 

}

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "NURBS_CANS");
  ros::NodeHandle nh;

  mp.numRowsDesired = 95;
  mp.numColsDesired = 95;
  mp.maxNanAllowed = 10;
  mp.removeNanBuffer = 3;
  mp.nCtrlDefault[0] = 30;
  mp.nCtrlDefault[1] = 30;

  // Create a ROS subscriber for the input point cloud
  ros::Subscriber sub = nh.subscribe ("cloud", 1, cloud_cb);

  
  // Spin
  ros::spin ();
}
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

#include <tf/transform_broadcaster.h>

// CANS includes
// #include "cans/mapping3D.h"
// #include "cans/object3D.h"
// #include "cans/splitSurface.h"

using namespace std;
using namespace PLib;


sensor_msgs::PointCloud2 loadScanFromFile(std::string& filestem, int& i){

  std::string filename;
  pcl::PCLPointCloud2 cloud_blob;
  pcl::PCDReader reader;

  // Get filename:
  if (i < 9){
    filename = filestem + "0" + static_cast<ostringstream*>( &(ostringstream() << (i+1)) )->str() + ".pcd";
    cout << "filename is : " << filename;
  }else{
    filename = filestem + static_cast<ostringstream*>( &(ostringstream() << (i+1)) )->str() + ".pcd";
    cout << "filename is : " << filename;
  }

  // Read Scan
  reader.read (filename, cloud_blob);

  // Convert to ROS data type 
  sensor_msgs::PointCloud2 output;
  pcl_conversions::moveFromPCL(cloud_blob, output);

  output.header.seq = i;
  output.header.stamp = ros::Time(0);

  return output;
}

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "CloudPublisher");
  ros::NodeHandle nh;
  ros::Publisher pub;

  // Create a ROS publisher for a point cloud
  pub = nh.advertise<sensor_msgs::PointCloud2> ("cloud", 1);

  std::string filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Data000";

  sensor_msgs::PointCloud2 cloud;  

  ros::Rate r(0.1); // 0.1 hz

  int i = 0;

  while (nh.ok() && i < 10){

    pub.publish(loadScanFromFile(filestem, i));

    ros::spinOnce();
    r.sleep();
    i++;
  }
  
  

}


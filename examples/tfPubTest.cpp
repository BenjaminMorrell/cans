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

using namespace std;

tf::Transform transform_current;
tf::Quaternion q;
Eigen::Array<float,6,10> state;
int i;

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

    state(0,i) = f2;
    state(1,i) = f3;
    state(2,i) = f4;
    state(3,i) = f5;
    state(4,i) = f6;
    state(5,i) = f7;
    
  }

  fclose(pFile);

  return state;
}

void 
cloud_cb (const sensor_msgs::PointCloud2ConstPtr& cloud_msg)
{
  
  transform_current.setOrigin( tf::Vector3(state(0,i), state(1,i), state(2,i)));
  q.setRPY(state(3,i), state(4,i), state(5,i));
  transform_current.setRotation(q);

  i++;
  cout << "i is " << i << endl;
}

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "tfPub");
  ros::NodeHandle nh;

  state = readPathTextFile("/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Path.txt");
  static tf::TransformBroadcaster br;
  cout << "State is:\n" << state << endl;
  

  // Create a ROS subscriber for the input point cloud
  ros::Subscriber sub = nh.subscribe ("cloud", 1, cloud_cb);

  ros::Rate r(10.0); //  hz

  i = 0;


  while (nh.ok()){

    
    br.sendTransform(tf::StampedTransform(transform_current, ros::Time::now(), "world", "body"));

    ros::spinOnce();
    r.sleep();
  }
  
  // Spin
  // ros::spin ();
}


  
  
  

  
  
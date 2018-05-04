#include <ros/ros.h>
#include <iostream>
#include <fstream>

// PCL specific includes
#include <sensor_msgs/PointCloud2.h>
#include <cans_msgs/Object3D.h>
#include <pcl_conversions/pcl_conversions.h>


#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include "tf/transform_datatypes.h"
#include "tf_conversions/tf_eigen.h"
#include "eigen_conversions/eigen_msg.h"

// CANS includes
#include "cans/nurbSLAM.h"


using namespace std;
using namespace PLib;

class nurbSLAMNode {
  public:

    NurbSLAM slam;

    pcl::PointCloud<pcl::PointNormal>::Ptr cloud;
    std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds;

    Eigen::Affine3d transformTruth;
    Eigen::Affine3f state;

    tf::StampedTransform transformTF;

    bool runFlag;
    bool useTruthTransform;
    bool bNewObjects;
    bool bNewState;
    bool bNewScanReceived;

    int scanNumber;

  public:
    //------------------------------------------
    // Constructor
    nurbSLAMNode(): runFlag(true), useTruthTransform(false), scanNumber(0), 
        bNewObjects(false), bNewState(false), cloud(new pcl::PointCloud<pcl::PointNormal>),
        bNewScanReceived(false)
    {
      state = Eigen::Affine3f::Identity();
      transformTruth = Eigen::Affine3d::Identity();

      slam.bShowAlignment = false;
      slam.alignmentOption = 0; // 0 - dense to dense, 1 - keypoints to dense, 2 - keypoints to keypoints
      slam.pclNormalRadiusSetting = 0.05;
      slam.pclFeatureRadiusSetting = 0.1;
      slam.nSurfPointsFactor = 5.0; // default is 3.0
      slam.localisationOption = 2;// Option for localisation method (0 - PCL, 1 - RANSAC IA, 2 - Prerejective RANSAC)

      // Other settings
      slam.ransac_inlierMultiplier = 0.1; // for the RANSAC inlier threshold
      slam.modelResolutionKeypoints = 0.005; // For localisation

      transformTF.child_frame_id_ = "nurb_cam";
      transformTF.frame_id_ = "starting_cam";
    }

    // Process a point cloud message 
    void cloud_cb(const sensor_msgs::PointCloud2ConstPtr& cloud_msg){
      // Don't run if already running (runFlag = false)
      if (runFlag){
        // Set flag to false until finished
        runFlag = false;
        cout << "\nProcessing Scan " << scanNumber << "./n/n";

        // Empty cloud vector
        clouds.clear();

        if (useTruthTransform){
          tf::TransformListener listener;
          tf::StampedTransform transform;

          try{
            // listener.waitForTransform("body", "world", cloud_msg->header.stamp, ros::Duration(0.1));
            // listener.lookupTransform("body", "world",  
            //                           cloud_msg->header.stamp, transform);
            listener.waitForTransform("cam_optical", "world", ros::Time(0), ros::Duration(0.5));
            listener.lookupTransform("cam_optical", "world",  
                                      ros::Time(0), transform);                          
          }
          catch (tf::TransformException ex){
            ROS_ERROR("%s",ex.what());
            // ros::Duration(1.0).sleep();
          }
          // Fill transform
          tf::transformTFToEigen(transform,transformTruth);

          cout << "transform is " << transformTruth.matrix() << endl;
        }
        // Container for original & filtered data
        pcl::PCLPointCloud2* cloud_blob = new pcl::PCLPointCloud2; 
        
        // Convert to PCL data type
        pcl_conversions::toPCL(*cloud_msg, *cloud_blob);

        // Convert to Point Cloud <T>
        pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 

        // Add to point cloud vector (for later potentail use with multiple objects)
        clouds.push_back(cloud);
        
        std::cerr << "PointCloud Received: size: " << cloud->width * cloud->height << endl;
        
        // Set flag back to true - ready for next observation
        runFlag = true;

        // Set flag that a new scan has been received
        bNewScanReceived = true;
        
      }
    }

    //------------------------------------------
    // Process the point cloud on a timer
    void processPointCloud(const ros::TimerEvent&){
      // Process on a timer to process the scan
      if (!bNewScanReceived){
        cout << "No new scans, not doing anything in processPointCloud" << endl;
        return;
      }

      // Set flag to false - are processing this scan
      bNewScanReceived = false;
      
      
      if (useTruthTransform && slam.isMappingModeActive()){
        // Set the state
        slam.setState(transformTruth.cast<float>()); 
      }

      // Process scans
      slam.processScans(clouds);
      
      // Get the state
      state = slam.getState();

      // Set flags that there has been an update
      bNewState = true;

      if (slam.bMapUpdatedFromScan){
        // If there have been updates
        bNewObjects = true; 
      }

      
      std::string filename = "testNURBS_Unreal_" + static_cast<ostringstream*>( &(ostringstream() << (scanNumber)) )->str() + ".wrl";
      slam.mp.objectMap[slam.mp.objectMap.size()-1].writeVRML(filename.c_str(),Color(255,100,255),50,80); 

      cout << "\nFinished Processing Scan " << scanNumber << ".\n\n";
      scanNumber ++;
    }

    //------------------------------------------
    // Fill object message for a give ID to send as a ROS message
    cans_msgs::Object3D fillObject3DMessageForID(int objID){

      // Initialise message
      cans_msgs::Object3D msg;

      if (objID >= slam.mp.objectMap.size()){
        cout << "Error, objID is larger than map size. Returning empty message" << endl;
        msg.ID = -1;
        return msg;
      }
      
      // Fill message
      slam.mp.fillObject3DMessage(objID, msg);

      // Set flag back to false - objects are no longer new
      bNewObjects = false;

      return msg;
    }

    //------------------------------------------
    // Fill transform to then send
    tf::StampedTransform getCurrentTransform(){
      // transformTF is a member variable with frames already defined

      // Convert from state to transformTF
      tf::transformEigenToTF(state.cast<double>(),transformTF);

      // Add timestamp
      transformTF.stamp_ = ros::Time::now();

      // Set flag back to false - state is no longer new
      bNewState = false;
      
      return transformTF;
    }

    void setSLAMParameters(ros::NodeHandle& nh){

      cout << "Inside set parameters" << endl;
      // OPTIONS
      nh.param("alignmentOption", slam.alignmentOption, slam.alignmentOption);
      nh.param("bShowAlignment", slam.bShowAlignment, slam.bShowAlignment);
      nh.param("localisationOption", slam.localisationOption, slam.localisationOption);
      nh.param("keypointOption", slam.keypointOption, slam.keypointOption);

      // Localisation
      nh.param("/keypoints/modelResolution", slam.modelResolutionKeypoints, slam.modelResolutionKeypoints);
      nh.param("/keypoints/minNeighbours", slam.minNeighboursKeypoints, slam.minNeighboursKeypoints);

      nh.param("pclNormalRadiusSetting", slam.pclNormalRadiusSetting, slam.pclNormalRadiusSetting);
      nh.param("pclFeatureRadiusSetting", slam.pclFeatureRadiusSetting, slam.pclFeatureRadiusSetting);

      nh.param("/ransac/inlierMultiplier", slam.ransac_inlierMultiplier, slam.ransac_inlierMultiplier);
      nh.param("/ransac/maximumIterations", slam.ransac_maximumIterations, slam.ransac_maximumIterations);
      nh.param("/ransac/numberOfSamples", slam.ransac_numberOfSamples, slam.ransac_numberOfSamples);
      nh.param("/ransac/correspondenceRandomness", slam.ransac_correspondenceRandomness, slam.ransac_correspondenceRandomness);
      nh.param("/ransac/similarityThreshold", slam.ransac_similarityThreshold, slam.ransac_similarityThreshold);
      nh.param("/ransac/inlierFraction", slam.ransac_inlierFraction, slam.ransac_inlierFraction);

      nh.param("validInlierThreshold", slam.validInlierTheshold, slam.validInlierTheshold);
      nh.param("nSurfPointsFactor", slam.nSurfPointsFactor, slam.nSurfPointsFactor);

      // Mapping
      nh.param("/meshing/numRowsDesired", slam.mp.numRowsDesired, slam.mp.numRowsDesired);
      nh.param("/meshing/numColsDesired", slam.mp.numColsDesired, slam.mp.numColsDesired);
      nh.param("/meshing/maxNanAllowed", slam.mp.maxNanAllowed, slam.mp.maxNanAllowed);
      nh.param("/meshing/removeNanBuffer", slam.mp.removeNanBuffer, slam.mp.removeNanBuffer);
      nh.param("/meshing/newRowColBuffer", slam.mp.newRowColBuffer, slam.mp.newRowColBuffer);
      nh.param("/meshing/useNonRectData", slam.mp.useNonRectData, slam.mp.useNonRectData);
      nh.param("/meshing/nCtrlDefaultS", slam.mp.nCtrlDefault[0], slam.mp.nCtrlDefault[0]); // This main fail...
      nh.param("/meshing/nCtrlDefaultT", slam.mp.nCtrlDefault[1], slam.mp.nCtrlDefault[1]);

      cout << "Finished setting SLAM parameters" << endl;

    }
};

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "NURBSLAM");
  ros::NodeHandle nh;

  // Initialise the class
  nurbSLAMNode nurbnode;
  nurbnode.setSLAMParameters(nh);

  if (argc > 1){
    int slamMode = atoi(argv[1]);
    switch (slamMode){
      case 0:
        cout << "\n\n\t\tACTIVATING SLAM MODE (default)\n\n" << endl;
        break;
      case 1:
        nurbnode.slam.activateMappingMode();
        cout << "\n\n\t\tACTIVATING PURE MAPPING MODE\n\n" << endl;
        break;
      case 2:
        nurbnode.slam.activateLocalisationMode();
        cout << "\n\n\t\tACTIVATING PURE LOCALISATION MODE\n\n" << endl;
        // load object
        try{
          throw 1;
          // nurbnode.slam.loadObjectIntoMap((outFilestem + "_0_final_obj.obj").c_str());
          // TODO may need to update this to load multiple objects
        }catch(...){
          cout << "\n\nLocalisation not yet implemented... No object found to use for localisation. Exiting." << endl;
          return -1;
        }
        break;
    }
  }

  // Create a ROS subscriber for the input point cloud
  ros::Subscriber sub = nh.subscribe ("/camera/points2", 1, &nurbSLAMNode::cloud_cb, &nurbnode);

  // Timer to process the cloud
  ros::Timer timer = nh.createTimer(ros::Duration(10.0), &nurbSLAMNode::processPointCloud, &nurbnode);
  
  // Publisher for the NURBS objects
  ros::Publisher mapPub = nh.advertise<cans_msgs::Object3D>("object",1,false);

  tf::TransformBroadcaster tf_br;

  // Spin
  // ros::spin ();
  
  ros::Rate r(2);// Adjust this 
  while (nh.ok()){
    if (nurbnode.bNewObjects){
      // If the objects have been updated
      for (int i = 0; i < nurbnode.slam.mp.objectMap.size(); i++){
        // for each object - fill a message and publish it
        mapPub.publish(nurbnode.fillObject3DMessageForID(i));
      }
    }
    if (true || nurbnode.bNewState){ // Or just publish on every loop
      tf_br.sendTransform(nurbnode.getCurrentTransform());
    }
    
    ros::spinOnce();
    r.sleep();
  }
}
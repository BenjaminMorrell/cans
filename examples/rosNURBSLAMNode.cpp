#include <ros/ros.h>
#include <ros/package.h>
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

    char* stateFilename;

    int scanNumber;

    
  public:
    //------------------------------------------
    // Constructor
    nurbSLAMNode(): runFlag(true), useTruthTransform(false), scanNumber(0), 
        bNewObjects(false), bNewState(false), cloud(new pcl::PointCloud<pcl::PointNormal>),
        bNewScanReceived(false), stateFilename("/home/bjm/SpaceCRAFT/unrealDataTrack.txt")
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

      // Clear file
      ofstream myfile;
      myfile.open (stateFilename, std::ofstream::out | std::ofstream::trunc); // open and close to clear.
      myfile.close();

    }

    ~nurbSLAMNode(){
      ;
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

        // pcl::PCDWriter writer;
        // writer.write<pcl::PointNormal> ("/home/amme2/Development/voxblox_ws/unreal_scan.pcd", *cloud, false);
        
      }
    }

    void updateStateFile(){

      Eigen::Array<float, 1, 6> stateData(1,6);

      // Update State for tracking
      stateData[0] = (state.matrix()(0,3));
      stateData[1] = (state.matrix()(1,3));
      stateData[2] = (state.matrix()(2,3));
      
      Eigen::Vector3f rpy; // init vector to store output state

      rpy = state.rotation().eulerAngles(2, 1, 0); // Check ordering

      stateData[3] = (rpy(2));
      stateData[4] = (rpy(1));
      stateData[5] = (rpy(0));
      

      // Open file
      ofstream myfile;
      myfile.open (stateFilename, std::ios_base::app); // Append
      myfile << stateData << "\n";
      myfile.close();
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

      // /home/amme2/Development/voxblox_ws/testNURBS_Unreal_6.wrl
      // TODO - wrap this in an if statement
      // for (int i=0; i < slam.mp.objectMap.size(); i++){
      //   std::string filename = "/home/amme2/Development/Results/testNURBS_Unreal_" + static_cast<ostringstream*>( &(ostringstream() << (i)) )->str() + ".wrl";
      //   slam.mp.objectMap[i].writeVRML(filename.c_str(),Color(255,100,255),50,80); 
      //   filename = "/home/amme2/Development/Results/testNURBS_Unreal_" + static_cast<ostringstream*>( &(ostringstream() << (i)) )->str() + ".pcd";
      //   slam.mp.writeObjectPCDFile(filename.c_str(), i, 125, 125);
      // }
      
      cout << "\nFinished Processing Scan " << scanNumber << ".\n\n";
      scanNumber ++;

      // Store state
      updateStateFile();
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

      // Set flag back to false - objects are no longer new (will still loop through all in the map)
      // bNewObjects = false;

      return msg;
    }

    //------------------------------------------
    // Fill object point cloud message for a give ID to send as a ROS message
    sensor_msgs::PointCloud2 fillPointCloudMessageForID(int objID){

      // Initialise message
      sensor_msgs::PointCloud2 msg;

      if (objID >= slam.mp.objectMap.size()){
        cout << "Error, objID is larger than map size. Returning empty message" << endl;
        return msg;
      }

      // Number of data points
      int ms = 30;
      int mt = 50;
      // Init point cloud
      pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>(mt, ms, pcl::PointNormal()));

      // Get point cloud from object
      slam.mp.pointCloudFromObject3D(objID, ms, mt, cloud);
      
      // Convert to ROS message (http://docs.ros.org/hydro/api/pcl_conversions/html/namespacepcl.html#af2c39730f92ade1603c55d45265e386d)
      pcl::toROSMsg(*cloud, msg);

      // Set frame ID
      msg.header.frame_id = "local_origin";
      msg.header.seq = objID;
      msg.header.stamp = ros::Time::now();

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
      nh.param("bRejectNonOverlappingInAlign", slam.bRejectNonOverlappingInAlign, slam.bRejectNonOverlappingInAlign);

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

      nh.param("maxDistanceOverlap", slam.maxDistanceOverlap, slam.maxDistanceOverlap);

      // Mapping
      nh.param("/meshing/numRowsDesired", slam.mp.numRowsDesired, slam.mp.numRowsDesired);
      nh.param("/meshing/numColsDesired", slam.mp.numColsDesired, slam.mp.numColsDesired);
      nh.param("/meshing/maxNanAllowed", slam.mp.maxNanAllowed, slam.mp.maxNanAllowed);
      nh.param("/meshing/removeNanBuffer", slam.mp.removeNanBuffer, slam.mp.removeNanBuffer);
      nh.param("/meshing/newRowColBuffer", slam.mp.newRowColBuffer, slam.mp.newRowColBuffer);
      nh.param("/mapping/useNonRectData", slam.mp.useNonRectData, slam.mp.useNonRectData);
      nh.param("/mapping/nCtrlDefaultS", slam.mp.nCtrlDefault[0], slam.mp.nCtrlDefault[0]); 
      nh.param("/mapping/nCtrlDefaultT", slam.mp.nCtrlDefault[1], slam.mp.nCtrlDefault[1]);

      cout << "nCtrlDefaultS is " << slam.mp.nCtrlDefault[0] << endl;
      cout << "nCtrlDefaultT is " << slam.mp.nCtrlDefault[1] << endl;

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

  bool bPublishNurbsPointCloud = false;
  nh.param("bPublishNurbsPointCloud", bPublishNurbsPointCloud, bPublishNurbsPointCloud);

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
          std::string filename = ros::package::getPath("cans")+"/data/blob_0_final_obj.obj";
          nurbnode.slam.loadObjectIntoMap(filename.c_str());
          filename = ros::package::getPath("cans")+"/data/blob_1_final_obj.obj";
          nurbnode.slam.loadObjectIntoMap(filename.c_str());
          filename = ros::package::getPath("cans")+"/data/blob_2_final_obj.obj";
          nurbnode.slam.loadObjectIntoMap(filename.c_str());
          filename = ros::package::getPath("cans")+"/data/blob_3_final_obj.obj";

          nurbnode.bNewObjects = true;
          // TODO may need to update this to load multiple objects
        }catch(...){
          cout << "\n\nFiles not present to do localistion" << endl;//Localisation not yet implemented... No object found to use for localisation. Exiting." << endl;
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
  ros::Publisher mapPCPub = nh.advertise<sensor_msgs::PointCloud2>("object_point_cloud",1,false);

  tf::TransformBroadcaster tf_br;

  // Spin
  // ros::spin ();
  
  ros::Rate r(1);// Adjust this 
  int counter = 0;
  while (nh.ok()){
    if (counter%5 == 0){//nurbnode.bNewObjects){
      // If the objects have been updated
      for (int i = 0; i < nurbnode.slam.mp.objectMap.size(); i++){
        // for each object - fill a message and publish it
        cout << "Publishing message for object with ID: " << i << endl;
        mapPub.publish(nurbnode.fillObject3DMessageForID(i));

        if (bPublishNurbsPointCloud){
          mapPCPub.publish(nurbnode.fillPointCloudMessageForID(i));
        }
      }
      // nurbnode.bNewObjects = false;
    }

    if (true || nurbnode.bNewState){ // Or just publish on every loop
      tf_br.sendTransform(nurbnode.getCurrentTransform());
    }

    counter++;
    
    ros::spinOnce();
    r.sleep();
  }
}
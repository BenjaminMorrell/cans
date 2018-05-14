#include <ros/ros.h>
#include <ros/package.h>
#include <iostream>
#include <fstream>

//PCL includes
#include <Eigen/Core>
// #include <pcl/point_types.h>
// #include <pcl/point_cloud.h>
// #include <pcl/common/time.h>

// #include <pcl/features/normal_3d_omp.h>
// #include <pcl/features/normal_3d.h>
// #include <pcl/features/fpfh_omp.h>

// #include <pcl/common/transforms.h>

// #include <pcl/io/pcd_io.h>
// #include <pcl/registration/icp.h>
// #include <pcl/registration/ia_ransac.h>
// #include <pcl/registration/sample_consensus_prerejective.h>
// #include <pcl/registration/correspondence_estimation_normal_shooting.h>
// #include <pcl/registration/correspondence_rejection_distance.h>
// #include <pcl/segmentation/sac_segmentation.h>
// #include <pcl/visualization/pcl_visualizer.h>

// #include <pcl/visualization/cloud_viewer.h>

// #include "cans/mapping3D.h"
// #include "cans/object3D.h"
#include "cans/nurbSLAM.h"

using namespace std;


Eigen::Array<float,6,Eigen::Dynamic> readPathTextFile(const char * filename, int nData){

  std::ifstream file;
  std::string line;

  Eigen::Array<float,6,Eigen::Dynamic> state(6,nData);

  float f1, f2, f3, f4, f5, f6, f7;
  char* s1[10], s2[7], s3[6], s4[2], s5[2], s6[2], s7[3], s8[3], s9[3];

  FILE * pFile;

  pFile = fopen (filename, "r");

  for (int i = 0; i < nData; i++){
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

void setSLAMParameters(NurbSLAM& slam, ros::NodeHandle nh){

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

  nh.param("mapCountThreshold", slam.mapCountThreshold, slam.mapCountThreshold);
  nh.param("mapExtendThreshold", slam.mapExtendThreshold, slam.mapExtendThreshold);
  

  

  // Mapping
  nh.param("/meshing/numRowsDesired", slam.mp.numRowsDesired, slam.mp.numRowsDesired);
  nh.param("/meshing/numColsDesired", slam.mp.numColsDesired, slam.mp.numColsDesired);
  nh.param("/meshing/maxNanAllowed", slam.mp.maxNanAllowed, slam.mp.maxNanAllowed);
  nh.param("/meshing/removeNanBuffer", slam.mp.removeNanBuffer, slam.mp.removeNanBuffer);
  nh.param("/meshing/newRowColBuffer", slam.mp.newRowColBuffer, slam.mp.newRowColBuffer);
  nh.param("/meshing/bFilterZ", slam.mp.bFilterZ, slam.mp.bFilterZ);
  nh.param("/meshing/nPointsZLim", slam.mp.nPointsZLim, slam.mp.nPointsZLim);
  nh.param("/meshing/bNegateZ", slam.mp.bNegateZ, slam.mp.bNegateZ);

  nh.param("/mapping/useNonRectData", slam.mp.useNonRectData, slam.mp.useNonRectData);
  nh.param("/mapping/nCtrlDefaultS", slam.mp.nCtrlDefault[0], slam.mp.nCtrlDefault[0]); 
  nh.param("/mapping/nCtrlDefaultT", slam.mp.nCtrlDefault[1], slam.mp.nCtrlDefault[1]);

  nh.param("/mapping/bUseFullAlignmentTransformInUpdate", slam.bUseFullAlignmentTransformInUpdate, slam.bUseFullAlignmentTransformInUpdate);
  nh.param("/mapping/bUseOldStateForNewObjects", slam.bUseOldStateForNewObjects, slam.bUseOldStateForNewObjects);
  

  cout << "nCtrlDefaultS is " << slam.mp.nCtrlDefault[0] << endl;
  cout << "nCtrlDefaultT is " << slam.mp.nCtrlDefault[1] << endl;

  // SLAM EKF
  nh.param("/ekf/pNoisePos", slam.pNoisePos, slam.pNoisePos);
  nh.param("/ekf/pNoiseVel", slam.pNoiseVel, slam.pNoiseVel);
  nh.param("/ekf/pNoiseAccel", slam.pNoiseAccel, slam.pNoiseAccel);
  nh.param("/ekf/pNoiseAng", slam.pNoiseAng, slam.pNoiseAng);
  nh.param("/ekf/pNoiseMultiplier", slam.pNoiseMultiplier, slam.pNoiseMultiplier);
  nh.param("/ekf/qNoiseMultiplier", slam.qNoiseMultiplier, slam.qNoiseMultiplier);
  
  nh.param("/ekf/noiseObsBasePos", slam.noiseObsBasePos, slam.noiseObsBasePos);
  nh.param("/ekf/noiseObsMultPos", slam.noiseObsMultPos, slam.noiseObsMultPos);
  nh.param("/ekf/noiseObsMultPosErr", slam.noiseObsMultPosErr, slam.noiseObsMultPosErr);
  nh.param("/ekf/noiseObsBaseAng", slam.noiseObsBaseAng, slam.noiseObsBaseAng);
  nh.param("/ekf/noiseObsMultAng", slam.noiseObsMultAng, slam.noiseObsMultAng);
  nh.param("/ekf/noiseObsMultAngErr", slam.noiseObsMultAngErr, slam.noiseObsMultAngErr);
  nh.param("/ekf/rMatMultiplier", slam.rMatMultiplier, slam.rMatMultiplier);
  nh.param("/ekf/bKeepPConstant", slam.bKeepPConstant, slam.bKeepPConstant);

  nh.param("/ekf/processModel", slam.processModel, slam.processModel);

  nh.param("/ekf/rejectCriteriaAng", slam.rejectCriteria[0], slam.rejectCriteria[0]);
  nh.param("/ekf/rejectCriteriaLin", slam.rejectCriteria[2], slam.rejectCriteria[2]);
  nh.param("/ekf/rejectCriteriaInlier", slam.rejectCriteria[4], slam.rejectCriteria[4]);
  nh.param("/ekf/rejectCriteriaNumberP", slam.rejectCriteria[5], slam.rejectCriteria[5]);
  slam.rejectCriteria[3] = 2.0*slam.rejectCriteria[2];
  slam.rejectCriteria[1] = 2.0*slam.rejectCriteria[0];
  


  cout << "Finished setting SLAM parameters" << endl;

}

void runSLAM(int argc,char ** argv, ros::NodeHandle nh){
  // Init
  int numberOfScans;
  int dataSet;
  int scanSteps;
  int nData;
  int nurbSLAMMode;

  if (argc < 3){
    cout << "Error: need to use at least 2 arugments: int dataset, int numberOfScans (optional) mode {0, SLAM (default), 1, mapping, 2 localisation}" << endl;
    return;
  }else{
    dataSet = atoi(argv[1]);
    numberOfScans = atoi(argv[2]);
    if (argc > 3){
      nurbSLAMMode = atoi(argv[3]);
    }
    
  }

  // Initialise
  std::string filename;
  std::string outFilename;
  std::string filestem;
  std::string outFilestem;
  std::string pathFilename;

  Eigen::Affine3f transform = Eigen::Affine3f::Identity();
  Eigen::Affine3f alignmentTransform = Eigen::Affine3f::Identity();

  // Points to generate on the surface

  // Load data for given test case:
  cout << "In blender sequence test " << endl;

  int machine = 0;
  nh.param("machine", machine, machine);

  switch (dataSet){
    case 0: 
      // BLOB
      cout << "Running Blob Dataset\n";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/blob/blob_scan_res_new_extend_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Path.txt";
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 10;
      break;
    case 1:
      // Longer Blob
      cout << "Running Long Blob Dataset\n";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/BlobLong/BlobScan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/BlobLong/blob_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/BlobLong/BlobScan_path.txt";
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 100;
      break;
    case 2:
      // Cube
      cout << "Running Cube Dataset\n";
      // filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan02/BlockScan_data00";
      // outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/Cube/cube_scan_res_";
      // pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan02/BlockScan_path.txt";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/CUBE/oneObj_scan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/Cube/cube_one_obj_scan_res_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/CUBE/oneObj_scan_path.txt";
      // Also CUBE_DIAG, CUBE_VERT, 
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 10;
      break;
    case 3:
      // Cube
      cout << "\n\t\tRunning Cube VERT Dataset\n\n";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/CUBE_VERT/oneObj_scan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/Cube/cube_vert_one_obj_scan_res_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/CUBE_VERT/oneObj_scan_path.txt";
      // Also CUBE_DIAG, CUBE_VERT, 
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 10;
      break;
    case 4:
      // Sphere
      cout << "Running Sphere Dataset\n";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/SPHERE/oneObj_scan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/Sphere/sphere_scan_res_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/SPHERE/oneObj_scan_path.txt";
      // ALso SPHERE_DIAG
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 10;
      break;
    case 5:
      // Sphere
      cout << "Running Sphere Dataset\n";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/SPHERE_DIAG/oneObj_scan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/Sphere/sphere_diag_scan_res_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/SPHERE_DIAG/oneObj_scan_path.txt";
      // ALso SPHERE_DIAG
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 10;
      break;
    case 6:
      // Longer Blob 2 - no problematic square parts
      if (machine == 0){
        filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Blob2/BlobScan_data00";
        outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/Blob2/blob_";
        pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Blob2/BlobScan_path.txt";
      }else if (machine == 1){
        filestem = "/home/amme2/Development/Data/Blob2/BlobScan_data00";
        outFilestem = "/home/amme2/Development/Results/Blensor/Blob2/blob_";
        pathFilename = "/home/amme2/Development/Data/Blob2/BlobScan_path.txt";
      }
      cout << "Running Long Blob2 Dataset\n";
      
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 100;
      break;
    case 7:
      // Blob Lateral - NOTE that the path may be wrong - so only use for SLAM...
      cout << "Running Long BlobLat Datasetn\n";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/BlobLat/BlobScan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/BlobLat/blob_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/BlobLat/BlobScan_path.txt";
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 40;
      break;
    case 99:
      // NEW BLOB
      cout << "Running New Blob Dataset\n";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/newBlob/newblobScans00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/newBlob/newblob_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/newBlob/newblobTrack.txt";
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 100;
      break;
  }

  // Init SLAM class
  NurbSLAM slam;

  slam.bShowAlignment = true;
  
  // slam.alignmentOption = 0; // 0 - dense to dense, 1 - keypoints to dense, 2 - keypoints to keypoints
  ros::param::get("alignmentOption", slam.alignmentOption);
  cout << "SLAM alignment option is " << slam.alignmentOption << endl;

  slam.pclNormalRadiusSetting = 0.05;
  slam.pclFeatureRadiusSetting = 0.1;
  slam.nSurfPointsFactor = 5.0; // default is 3.0
  slam.localisationOption = 2;// Option for localisation method (0 - PCL, 1 - RANSAC IA, 2 - Prerejective RANSAC)
  if (argc > 3){
    slam.localisationOption = atoi(argv[3]);
  }

  // Other settings
  // slam.ransac_inlierMultiplier = 0.1; // for the RANSAC inlier threshold
  // slam.modelResolutionKeypoints = 0.005; // For localisation

  setSLAMParameters(slam, nh);
  slam.setInitEKFStates();

  // Mapping settings
  // slam.mp.numRowsDesired = 95;
  // slam.mp.numColsDesired = 95;

  switch (nurbSLAMMode){
    case 0:
      cout << "\n\n\t\tACTIVATING SLAM MODE\n\n" << endl;
      break;
    case 1:
      slam.activateMappingMode();
      cout << "\n\n\t\tACTIVATING PURE MAPPING MODE\n\n" << endl;
      break;
    case 2:
      slam.activateLocalisationMode();
      cout << "\n\n\t\tACTIVATING PURE LOCALISATION MODE\n\n" << endl;

      break;
  }


  // Load true state
  Eigen::Array<float,6,Eigen::Dynamic> state = readPathTextFile(pathFilename.c_str(),nData);
  cout << "State is:\n" << state << endl;
  
  // Translation
  transform(0,3) = state(0,0);
  transform(1,3) = state(1,0);
  transform(2,3) = state(2,0);
  
    
  // 3, 2, 1 Euler transformation
  transform.rotate (Eigen::AngleAxisf(state(5,0),Eigen::Vector3f::UnitZ()));
  transform.rotate (Eigen::AngleAxisf(state(4,0),Eigen::Vector3f::UnitY()));
  transform.rotate (Eigen::AngleAxisf(state(3,0),Eigen::Vector3f::UnitX()));

  if (nurbSLAMMode==2){
    transform.setIdentity();
    transform(0,3) = 1.03987;
    transform(1,3) = 3.49792;
    transform(2,3) = 0.178754;

    transform.rotate (Eigen::AngleAxisf(3.0279,Eigen::Vector3f::UnitZ()));
    transform.rotate (Eigen::AngleAxisf(0.0100855,Eigen::Vector3f::UnitY()));
    transform.rotate (Eigen::AngleAxisf(1.67344,Eigen::Vector3f::UnitX()));
  }
  

  cout << "transform is " << transform.matrix() << endl;

  // Reset state to zero
  state.setZero(state.rows(),state.cols());

  // INITIALISE STATE
  slam.setState(transform);

  if (nurbSLAMMode==2){
    // load object
    try{
      slam.loadObjectIntoMap((outFilestem + "1_final_obj_saved.obj").c_str());
      // TODO may need to update this to load multiple objects
    }catch(...){
      cout << "\n\nNo object found to use for localisation. Exiting." << endl;
      return;
    }
  }

  // TIMESTEP
  float timestep = 0.1; 
  nh.param("/ekf/timestep", timestep, timestep);

  Eigen::Vector3f rpy; // init vector to store output state

  // Init pcl reader
  pcl::PCDReader reader;

  // Init clouds
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds;

  // init timers
  std::chrono::high_resolution_clock::time_point startTime;
  std::chrono::high_resolution_clock::time_point endTime;
  std::chrono::duration<double, std::milli> runtimeDuration;
  std::vector<double> processTimesVec(6);
  std::string timingFilename = ros::package::getPath("cans")+"/data/nurbsTimes.txt";
 
  // Clear file
  ofstream myfile;
  myfile.open (timingFilename.c_str(), std::ofstream::out | std::ofstream::trunc); // open and close to clear.
  myfile.close();
  ofstream myfile2;
  myfile2.open ((outFilestem + "state_track_slam.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
  myfile.close();
  
  int start_i = 0;

  if (nurbSLAMMode==2){
    start_i = 50;
  }

  // Run sequence of scans 
  for (int i = start_i; i < numberOfScans; i += scanSteps){
    cout << "Processing Scan " << i << endl;

    // if (i == 77){
    //   slam.bShowAlignment = true;
    // }

    // if (i == 60){
    //   slam.bShowAlignment = false;
    // }

    // Get filename:
    if (i < 9){
      filename = filestem + "00" + static_cast<ostringstream*>( &(ostringstream() << (i+1)) )->str() + ".pcd";
      cout << "filename is : " << filename;
    }else if (i < 99){
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
    cout << "converted PC" << endl;

    // Put in vector (if we had multiple objects)
    clouds.push_back(cloud);

    // Start timer
    startTime = std::chrono::high_resolution_clock::now(); 

    // Process scans
    slam.processScans(clouds, timestep);

    // End time and duraction
    endTime = std::chrono::high_resolution_clock::now();
    runtimeDuration = endTime - startTime;
    for (int i = 0; i < 5; i++){processTimesVec[i] = slam.processTimes[i];}
    processTimesVec[5] = runtimeDuration.count();

    // Get the state
    transform = slam.getState();

    // Update State for tracking
    state(0,i) = transform.matrix()(0,3);
    state(1,i) = transform.matrix()(1,3);
    state(2,i) = transform.matrix()(2,3);

    rpy = transform.rotation().eulerAngles(2, 1, 0); // Check ordering

    state(3,i) = rpy(2);
    state(4,i) = rpy(1);
    state(5,i) = rpy(0);

    // Empty cloud vector
    clouds.clear();

    // Write timing information
    myfile.open (timingFilename.c_str(), std::ios_base::app); // Append
    myfile << processTimesVec[0];
    for (int k=1; k < 6; k++){myfile << ", " << processTimesVec[k];}
    myfile << "\n";
    myfile.close();

    // Save final objects
    for (int j = 0; j < slam.mp.objectMap.size(); j++){
      
      // Write result to pcd
      filename = outFilestem + static_cast<ostringstream*>( &(ostringstream() << (j)) )->str() + "_end_result.pcd";
      slam.mp.writeObjectPCDFile(filename.c_str(), j, 125, 125);


      // Write final object to file
      filename = outFilestem + static_cast<ostringstream*>( &(ostringstream() << (j)) )->str() + "_final_obj.obj";
      slam.mp.objectMap[j].write(filename.c_str());

      // Write VRML
      filename = outFilestem + static_cast<ostringstream*>( &(ostringstream() << (j)) )->str() + "_savedObject.wrl";
      slam.mp.objectMap[j].writeVRML(filename.c_str(),Color(255,100,255),50,80); 
    }

    // Open file
    myfile.open ((outFilestem + "state_track_slam.txt").c_str(), std::ios_base::app); // Append
    myfile << state(0,i);
    for (int k=1; k < 6; k++){myfile << ", " << state(k,i);}
    myfile << "\n";
    myfile.close();
  }

  // Print final state
  cout << "Final state set is: " << state << endl;

  // Write to file
  myfile.open ((outFilestem + "state_track_slam_final.txt").c_str());
  myfile << state;
  myfile.close();

  
  
}



int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "nurbsSLAM_Testing");
  ros::NodeHandle nh;

  runSLAM(argc, argv, nh);
  // Spin
  // ros::spin ();
}
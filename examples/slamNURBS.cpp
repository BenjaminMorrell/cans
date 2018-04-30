#include <ros/ros.h>
#include <iostream>

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


void runSLAM(int argc,char ** argv){
  // Init
  int numberOfScans;
  int dataSet;
  int scanSteps;
  int nData;

  if (argc < 3){
    cout << "Error: need to use at least 2 arugments: int dataset, int numberOfScans (optional) int localisationMethod" << endl;
    return;
  }else{
    dataSet = atoi(argv[1]);
    numberOfScans = atoi(argv[2]);
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

  switch (dataSet){
    case 0: 
      // BLOB
      cout << "Running Blob Dataset";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/blob/blob_scan_res_new_extend_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Path.txt";
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 10;
      break;
    case 1:
      // Longer Blob
      cout << "Running Long Blob Dataset";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/BlobLong/BlobScan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/BlobLong/blob_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/BlobLong/BlobScan_path.txt";
      scanSteps = 5;
      numberOfScans = numberOfScans*scanSteps;
      nData = 100;
      break;
    case 2:
      // Cube
      cout << "Running Cube Dataset";
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
      cout << "Running Sphere Dataset";
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
      cout << "Running Sphere Dataset";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/SPHERE_DIAG/oneObj_scan_data00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/Sphere/sphere_diag_scan_res_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/One_Object/SPHERE_DIAG/oneObj_scan_path.txt";
      // ALso SPHERE_DIAG
      scanSteps = 1;
      numberOfScans = numberOfScans*scanSteps;
      nData = 10;
      break;
    case 99:
      // NEW BLOB
      cout << "Running New Blob Dataset";
      filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/newBlob/newblobScans00";
      outFilestem = "/home/bjm/Dropbox/PhD_Code/Results/newBlob/newblob_";
      pathFilename = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/newBlob/newblobTrack.txt";
      scanSteps = 5;
      numberOfScans = numberOfScans*scanSteps;
      nData = 100;
      break;
  }

  // Init SLAM class
  NurbSLAM slam;

  slam.bShowAlignment = true;
  slam.pclRadiusSetting = 0.2;
  slam.nSurfPointsFactor = 3.0; // default is 3.0
  slam.localisationOption = 2;// Option for localisation method (0 - PCL, 1 - RANSAC IA, 2 - Prerejective RANSAC)
  if (argc > 3){
    slam.localisationOption = atoi(argv[3]);
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

  cout << "transform is " << transform.matrix() << endl;

  // INITIALISE STATE
  slam.initState(transform);

  Eigen::Vector3f rpy; // init vector to store output state

  // Init pcl reader
  pcl::PCDReader reader;

  // Init clouds
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds;

  // Run sequence of scans 
  for (int i = 0; i < numberOfScans; i += scanSteps){
    cout << "Processing Scan " << i << endl;

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

    // Process scans
    slam.processScans(clouds);

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
  }
}


int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "nurbsSLAM_Testing");
  ros::NodeHandle nh;

  runSLAM(argc, argv);
  // Spin
  // ros::spin ();
}
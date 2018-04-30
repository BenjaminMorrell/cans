#include <ros/ros.h>
#include <iostream>

//PCL includes
#include <Eigen/Core>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/time.h>

#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/fpfh_omp.h>

#include <pcl/common/transforms.h>

#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/ia_ransac.h>
#include <pcl/registration/sample_consensus_prerejective.h>
#include <pcl/registration/correspondence_estimation_normal_shooting.h>
#include <pcl/registration/correspondence_rejection_distance.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/visualization/pcl_visualizer.h>

#include "cans/mapping3D.h"
#include "cans/object3D.h"
#include "cans/splitSurface.h"

using namespace std;

typedef pcl::visualization::PointCloudColorHandlerCustom<pcl::PointNormal> ColorHandlerT; // Visualisation?


// template <typename PointSource, typename PointTarget, typename FeatureT>
// class SampleConsensusPrerejective_Exposed : public pcl::SampleConsensusPrerejective<PointSource, PointTarget, FeatureT> {
//   public:
//     pcl::CorrespondencesPtr getCorrespondencesPtr() {
//       std::cout << "Inside getCorrespondencesPtr()" << std::endl;
//       std::cout << "Correspondences size is: " << this->correspondences_->size() << std::endl;
//       pcl::Correspondence currentCorrespondence1 = (*this->correspondences_)[0];
//       std::cout << "Got correspondance 0" << std::endl;
//       for (int i = 0; i < this->correspondences_->size(); i++){
//         pcl::Correspondence currentCorrespondence = (*this->correspondences_)[i];
//         std::cout << "Source Point Index: " << currentCorrespondence.index_query << std::endl;
//         std::cout << "Target Point Match Index: " << currentCorrespondence.index_match << std::endl;
//         std::cout << "Distance between correspondences: " << currentCorrespondence.distance << std::endl;
//         std::cout << "Confidence in correspondence: " << currentCorrespondence.weight << std::endl;
//       }
//       std::cout << "Outside loop" << endl;
//       return this->correspondences_;
//     }
// };
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


void alignScans(int argc, char** argv){

  // Perturb cloud 2
  if (argc != 9 && argc != 3 ){
    pcl::console::print_error ("Syntax is: %s cloudIn perturb_x, perturb_y, perturb_z theta vx vy vz\nOr: %s cloud1 cloud1", argv[0]);
    return;
  }
  // atof(argv[1]); // Converts argument to a float
  // cout << "input arguments are: " << atof(argv[1]) << ", " << atof(argv[2]) << ", " << atof(argv[3]) << endl;

  // pcl::PointXYZ perturb(atof(argv[1]),atof(argv[2]),atof(argv[3]));

  // Load the Point Cloud
  // pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2_orig (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2 (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2_aligned (new pcl::PointCloud<pcl::PointNormal>);

  pcl::PointCloud<pcl::FPFHSignature33>::Ptr cloud_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr cloud2_features (new pcl::PointCloud<pcl::FPFHSignature33>);

  // Fill in the cloud data
  pcl::PCDReader reader;
//   reader.read ("/home/bjm/SpaceCRAFT/ros_ws/src/pcl_testing/data/table_scene_lms400.pcd", *cloud_blob);
  // reader.read ("/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Data00010.pcd", *cloud_blob);
  // reader.read ("/home/bjm/Dropbox/PhD_Code/NURBS_2018_01_02/Test_pcd_00001.pcd", *cloud_blob);
  reader.read (argv[1], *cloud_blob);
  if (argc == 3){
    reader.read (argv[2], *cloud_blob2);
  }
  

  // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
  if (argc > 3){
    pcl::fromPCLPointCloud2 (*cloud_blob, *cloud2_orig); 
  }else{
    pcl::fromPCLPointCloud2 (*cloud_blob2, *cloud2); 
  }

  std::cerr << "PointCloud dimensions, W: " << cloud->width << "\tH: " << cloud->height << "\t data points." << std::endl;

  // Transform if there are inputs to do so
  if (argc > 3){
    // Create transformation matrix
    Eigen::Affine3f transform = Eigen::Affine3f::Identity();
    // Translation
    transform.translation() << atof(argv[2]), atof(argv[3]), atof(argv[4]);

    // Rotation
    // Rot vec
    Eigen::Vector3f rotVec(atof(argv[6]),atof(argv[7]),atof(argv[8]));
    rotVec.normalize();

    transform.rotate (Eigen::AngleAxisf (atof(argv[5]), rotVec));
    // transform.rotate (Eigen::Quaternionf q(2, 0, 1, -3)); 

    cout << "transformation matrix is: \n" << transform.matrix() << endl;
    // Will transform then rotate. i.e. multiply the 4x4 by [x;y;z;1] - get transformed x, y, z + dx, dy, dz

    // Transform
    pcl::transformPointCloud(*cloud2_orig, *cloud2, transform);

    // // Linear Perturbation
    // for (int i = 0; i < cloud2->width; i++){
    //   for (int j = 0; j < cloud2->height; j++){
    //     // cloud2->at(j,i).getArray3fMap() = cloud2->at(j,i).getArray3fMap() + perturb.getArray3fMap()

    //     cloud2->at(j,i).x = cloud2->at(j,i).x + perturb.x;
    //     cloud2->at(j,i).y = cloud2->at(j,i).y + perturb.y;
    //     cloud2->at(j,i).z = cloud2->at(j,i).z + perturb.z;
    //   }
    // }
  }

  // Estimate normals
  pcl::NormalEstimationOMP<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(0.025);
  nest.setInputCloud(cloud);
  nest.compute(*cloud);
  nest.setInputCloud(cloud2);
  nest.compute(*cloud2);

  // Estimate features
  pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  fest.setRadiusSearch(0.075);
  fest.setInputCloud(cloud);
  fest.setInputNormals(cloud);
  fest.compute (*cloud_features);
  fest.setInputCloud(cloud2);
  fest.setInputNormals(cloud2);
  fest.compute (*cloud2_features);

  // Alignment
  // SampleConsensusPrerejective_Exposed<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  pcl::SampleConsensusPrerejective<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  align.setInputSource(cloud2);
  align.setSourceFeatures(cloud2_features);
  align.setInputTarget(cloud);
  align.setTargetFeatures(cloud_features);
  //Settings
  align.setMaximumIterations (50000); // Number of RANSAC iterations
  align.setNumberOfSamples (3); // Number of points to sample for generating/prerejecting a pose
  align.setCorrespondenceRandomness (3); // Number of nearest features to use
  align.setSimilarityThreshold (0.9f); // Polygonal edge length similarity threshold
  align.setMaxCorrespondenceDistance (2.5f * 0.01f);// Inlier threshold
  align.setInlierFraction (0.75f); // Required inlier fraction for accepting a pose hypothesis
  //Perform alignement
  {
    pcl::ScopeTime t("Alignment");
    align.align (*cloud2_aligned);
  }

  // pcl::CorrespondencesPtr match_indices = align.getCorrespondencesPtr();

  // pcl::Correspondence currentCorrespondence = (*match_indices)[0];
  // std::cout << "Source Point Index: " << currentCorrespondence.index_query << std::endl;
  // std::cout << "Target Point Match Index: " << currentCorrespondence.index_match << std::endl;
  // std::cout << "Distance between correspondences: " << currentCorrespondence.distance << std::endl;
  // std::cout << "Confidence in correspondence: " << currentCorrespondence.weight << std::endl;

  // cout << (*match_indices)[0].index_query << endl;

  // cout << "\n\nPoint on srf 1: " << (*match_indices)[5].index_query << ", matches to point on srf 2: " << (*match_indices)[5].index_match << ", with distance: " << (*match_indices)[5].distance << std::endl;

  if (align.hasConverged()){
    printf("\n");
    Eigen::Matrix4f transformation = align.getFinalTransformation();
    cout << "transformation matrix is: \n" << transformation << endl;
    // pcl::console::print_info ("    | %6.3f %6.3f %6.3f | \n", transformation (0,0), transformation (0,1), transformation (0,2));
    // pcl::console::print_info ("    | %6.3f %6.3f %6.3f | \n", transformation (1,0), transformation (1,1), transformation (1,2));
    // pcl::console::print_info ("    | %6.3f %6.3f %6.3f | \n", transformation (2,0), transformation (2,1), transformation (2,2));
    // pcl::console::print_info ("\n");
    // pcl::console::print_info ("t = < %0.3f, %0.3f, %0.3f >\n", transformation (0,3), transformation (1,3), transformation (2,3));
    pcl::console::print_info ("\n");
    pcl::console::print_info ("Inliers: %i/%i\n", align.getInliers ().size (), cloud2->size ());
      
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    visu.addPointCloud (cloud, ColorHandlerT (cloud, 0.0, 255.0, 0.0), "cloud1");
    visu.addPointCloud (cloud2_aligned, ColorHandlerT (cloud2_aligned, 0.0, 0.0, 255.0), "cloud2_aligned");
    visu.addPointCloud (cloud2, ColorHandlerT (cloud2, 255.0, 0.0, 0.0), "cloud2");
    visu.spin ();

  // // Write the aligned data to disk
  // pcl::PCDWriter writer;
  // writer.write<pcl::PointNormal> ("pcl_aligned.pcd", *cloud2_aligned, false);

  }else{
    pcl::console::print_error ("Alignment failed!\n");
  }

}


void computePointCorrespondences(int argc, char ** argv){
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2 (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2_aligned (new pcl::PointCloud<pcl::PointNormal>);


  // Fill in the cloud data
  pcl::PCDReader reader;
  reader.read (argv[1], *cloud_blob);
  reader.read (argv[2], *cloud_blob2);

  // reader.read ("/home/bjm/Dropbox/PhD_Code/NURBS_2018_01_02/Test_pcd_1_global.pcd", *cloud_blob);
  // reader.read ("/home/bjm/Dropbox/PhD_Code/NURBS_2018_01_02/Test_pcd_2_global.pcd", *cloud_blob2);


  // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
  pcl::fromPCLPointCloud2 (*cloud_blob2, *cloud2); 

  // Estimate normals
  pcl::NormalEstimationOMP<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(0.025);
  nest.setInputCloud(cloud);
  nest.compute(*cloud);
  nest.setInputCloud(cloud2);
  nest.compute(*cloud2);

  // Estimate correspondences
  pcl::registration::CorrespondenceEstimationNormalShooting<pcl::PointNormal, pcl::PointNormal, pcl::PointNormal> corrEst;
  corrEst.setInputSource (cloud2);// Gets a corresponding point for every point in the Source
  corrEst.setSourceNormals (cloud2);
  corrEst.setInputTarget (cloud); // Target for the source to find corresponding points in
  // Test the first 10 correspondences for each point in source, and return the best
  corrEst.setKSearch (10);
  
  pcl::CorrespondencesPtr all_correspondencesPtr (new pcl::Correspondences);
  pcl::Correspondences& all_correspondences = *all_correspondencesPtr;
  
  // Determine all correspondences
  corrEst.determineCorrespondences (*all_correspondencesPtr);
  // corrEst.determineCorrespondences (all_correspondences);

  cout << " Correspondences determined " << endl;


  cout << "Correspondence 1, index query: " << all_correspondences[0].index_query << endl;
  cout << "Correspondence 1, match index: " << all_correspondences[0].index_match << endl;
  cout << "Correspondence 1, distance: " << all_correspondences[0].distance << endl;
  cout << "Correspondence 1, weight: " << all_correspondences[0].weight << endl;

  cout << "Correspondence 190, index query: " << all_correspondences[190].index_query << endl;
  cout << "Correspondence 190, match index: " << all_correspondences[190].index_match << endl;
  cout << "Correspondence 190, distance: " << all_correspondences[190].distance << endl;
  cout << "Correspondence 190, weight: " << all_correspondences[190].weight << endl;

  // Correspondence rejection
  pcl::CorrespondencesPtr corr_filtPtr (new pcl::Correspondences);
  pcl::Correspondences& corr_filt = *corr_filtPtr;
  pcl::registration::CorrespondenceRejectorDistance corrRej;
  corrRej.setMaximumDistance( atof(argv[3])); 
  corrRej.setInputCorrespondences( all_correspondencesPtr);
  corrRej.getCorrespondences( *corr_filtPtr);

  cout << "Index " << corr_filt[190].index_query << " matches to index " << corr_filt[190].index_match << endl;
  cout << "match weight: " << corr_filt[190].weight << endl;

  // -1 indicates that the correspondence has been rejected by the filter

  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> newPoints(cloud2->width,cloud2->height);
  newPoints.setZero(cloud2->width,cloud2->height);

  int offset;

  for (int i = 0; i < cloud2->height; i++){
    for (int j = 0; j < cloud2->width; j++){
      offset = i*cloud2->width + j; 
      if (corr_filt[offset].index_match == -1){
        // Not a correspondence, hence is a new point
        newPoints(i,j) = true;
      }
      // else leave at zero
    }
  }

  cout << "Number of new points: " << newPoints.count() << endl;

  Eigen::Array<int, Eigen::Dynamic, 1> newPointsRows = newPoints.rowwise().count().cast<int>();
  Eigen::Array<int, 1, Eigen::Dynamic> newPointsCols = newPoints.colwise().count().cast<int>();

  cout << "New count rows:\n" << newPointsRows << endl;
  cout << "New count cols:\n" << newPointsCols << endl;


  


}


void getAlignedCorrespondences(int argc, char ** argv){
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PCLPointCloud2::Ptr cloud_blob2 (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2 (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2_aligned (new pcl::PointCloud<pcl::PointNormal>);

  pcl::PointCloud<pcl::FPFHSignature33>::Ptr cloud_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr cloud2_features (new pcl::PointCloud<pcl::FPFHSignature33>);


  // Fill in the cloud data
  pcl::PCDReader reader;
  reader.read (argv[1], *cloud_blob);
  reader.read (argv[2], *cloud_blob2);

  // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud); 
  pcl::fromPCLPointCloud2 (*cloud_blob2, *cloud2); 

  // Estimate normals
  pcl::NormalEstimationOMP<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(0.025);
  nest.setInputCloud(cloud);
  nest.compute(*cloud);
  nest.setInputCloud(cloud2);
  nest.compute(*cloud2);

  // Estimate features
  pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  fest.setRadiusSearch(0.075);
  fest.setInputCloud(cloud);
  fest.setInputNormals(cloud);
  fest.compute (*cloud_features);
  fest.setInputCloud(cloud2);
  fest.setInputNormals(cloud2);
  fest.compute (*cloud2_features);

  // Alignment
  // SampleConsensusPrerejective_Exposed<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  pcl::SampleConsensusPrerejective<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  align.setInputSource(cloud2);
  align.setSourceFeatures(cloud2_features);
  align.setInputTarget(cloud);
  align.setTargetFeatures(cloud_features);
  //Settings
  align.setMaximumIterations (50); // Number of RANSAC iterations
  align.setNumberOfSamples (3); // Number of points to sample for generating/prerejecting a pose
  align.setCorrespondenceRandomness (3); // Number of nearest features to use
  align.setSimilarityThreshold (0.9f); // Polygonal edge length similarity threshold
  align.setMaxCorrespondenceDistance (2.5f * 0.01f);// Inlier threshold
  align.setInlierFraction (0.25f); // Required inlier fraction for accepting a pose hypothesis
  //Perform alignement
  {
    pcl::ScopeTime t("Alignment");
    align.align (*cloud2_aligned);
  }

  cout << "Completed alignment" << endl;
  // Get correspondences
  pcl::registration::CorrespondenceEstimationNormalShooting<pcl::PointNormal, pcl::PointNormal, pcl::PointNormal> corrEst;
  corrEst.setInputSource (cloud2_aligned);// Gets a corresponding point for every point in the Source
  corrEst.setSourceNormals (cloud2_aligned);
  corrEst.setInputTarget (cloud); // Target for the source to find corresponding points in
  // Test the first 10 correspondences for each point in source, and return the best
  corrEst.setKSearch (10);
  
  pcl::CorrespondencesPtr all_correspondencesPtr (new pcl::Correspondences);
  pcl::Correspondences& all_correspondences = *all_correspondencesPtr;
  
  cout << "Initialised correspondences" << endl;

  // Determine all correspondences
  corrEst.determineCorrespondences (*all_correspondencesPtr);

  cout << " Correspondences determined " << endl;

  cout << "Correspondence 1, index query: " << all_correspondences[0].index_query << endl;
  cout << "Correspondence 1, match index: " << all_correspondences[0].index_match << endl;
  cout << "Correspondence 1, distance: " << all_correspondences[0].distance << endl;
  cout << "Correspondence 1, weight: " << all_correspondences[0].weight << endl;

  cout << "Correspondence 190, index query: " << all_correspondences[190].index_query << endl;
  cout << "Correspondence 190, match index: " << all_correspondences[190].index_match << endl;
  cout << "Correspondence 190, distance: " << all_correspondences[190].distance << endl;
  cout << "Correspondence 190, weight: " << all_correspondences[190].weight << endl;

  // Correspondence rejection
  pcl::CorrespondencesPtr corr_filtPtr (new pcl::Correspondences);
  pcl::Correspondences& corr_filt = *corr_filtPtr;
  pcl::registration::CorrespondenceRejectorDistance corrRej;
  corrRej.setMaximumDistance( 0.01f); 
  corrRej.setInputCorrespondences( all_correspondencesPtr);
  corrRej.getCorrespondences( *corr_filtPtr);

  cout << "Rejected correspondences" << endl;

  cout << "Index " << corr_filt[190].index_query << " matches to index " << corr_filt[190].index_match << endl;
  cout << "match weight: " << corr_filt[190].weight << "\nMatch distance: " << corr_filt[190].distance << endl;

}

Eigen::Matrix4f runAlignmentPreRejective(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloud2, int option, bool bShowAlignment = true){
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr cloud_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr cloud2_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud2_aligned (new pcl::PointCloud<pcl::PointNormal>);

  cout << "inside runAlignmentPreRejective" << endl;

  cloud->is_dense = false;
  cloud2->is_dense = false;

  std::cerr << "PointCloud1 dimensions, W: " << cloud->width << "\tH: " << cloud->height << "\t is dense? " << cloud->is_dense << std::endl;
  std::cerr << "PointCloud2 dimensions, W: " << cloud2->width << "\tH: " << cloud2->height << "\t is dense? " << cloud2->is_dense << std::endl;

  

  // bool useicp = false;
  // bool useiaransac = false;
  // bool useprerejsac = true;

  Eigen::Matrix4f transformOut;

  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  if (option == 0){
    // ICP
    cout << "Starting alignment" << endl;
    pcl::IterativeClosestPoint<pcl::PointNormal, pcl::PointNormal> icp;
    // Set the input source and target
    icp.setInputSource (cloud2);
    icp.setInputTarget (cloud);
    // Set the max correspondence distance to 5cm (e.g., correspondences with higher distances will be ignored)
    icp.setMaxCorrespondenceDistance (0.5);
    // Set the maximum number of iterations (criterion 1)
    icp.setMaximumIterations (500);
    // Set the transformation epsilon (criterion 2)
    icp.setTransformationEpsilon (1e-8);
    // Set the euclidean distance difference epsilon (criterion 3)
    icp.setEuclideanFitnessEpsilon (1);
    
    // Perform the alignment
    icp.align (*cloud2_aligned);
    
    // Obtain the transformation that aligned cloud_source to cloud_source_registered
    transformOut = icp.getFinalTransformation ();
  }else{
    
    // Estimate normals
    pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
    nest.setRadiusSearch(0.1);
    // nest.setKSearch(7);
    nest.setSearchMethod(search_method_);
    nest.setInputCloud(cloud);
    nest.compute(*cloud);
    cout << "computed normals for cloud 1" << endl;
    nest.setInputCloud(cloud2);
    nest.compute(*cloud2);
    cout << "computed normals for cloud 2" << endl;
    // pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
    // nest.setRadiusSearch(0.01);
    // nest.setKSearch(5);
    // nest.setSearchMethod(pcl::search::KdTree<pcl::PointNormal>)
    // nest.setInputSource(cloud);
    // nest.compute(*cloud);
    // cout << "computed normals for cloud 1" << endl;
    // nest.setInputSource(cloud2);
    // nest.compute(*cloud2);

    cout << "Normals have been estimated" << endl;

    // Estimate features
    pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
    fest.setRadiusSearch(0.1);
    // fest.setKSearch(7);
    fest.setInputCloud(cloud);
    fest.setSearchMethod(search_method_);
    fest.setInputNormals(cloud);
    fest.compute (*cloud_features);
    fest.setInputCloud(cloud2);
    fest.setInputNormals(cloud2);
    fest.compute (*cloud2_features);
    cout << "Features have been computed" << endl;
  

    //Perform alignment 
    if (option == 1){
      // RANSAC Initial Alignment
      pcl::SampleConsensusInitialAlignment<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
      align.setInputSource(cloud2);
      align.setSourceFeatures(cloud2_features);
      align.setInputTarget(cloud);
      align.setTargetFeatures(cloud_features);
      //Settings
      align.setMaximumIterations (500); // Number of RANSAC iterations (1000 is default)
      align.setNumberOfSamples (3); // Number of points to sample for generating/prerejecting a pose (3 is default)
      align.setCorrespondenceRandomness (20); // Number of nearest features to use (default is 10)
      align.setMaxCorrespondenceDistance (0.1f);// Inlier threshold

      align.align(*cloud2_aligned);

      transformOut = align.getFinalTransformation ();

      if (align.hasConverged()){
        cout << " Successfully converged" << endl;
      }else {
        cout << "\n\tAlignment FAILED to converge\n" << endl;
      }
    }else if (option == 2){
      // Prerejective RANSAC
      // SampleConsensusPrerejective_Exposed<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
      pcl::SampleConsensusPrerejective<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
      align.setInputSource(cloud2);
      align.setSourceFeatures(cloud2_features);
      align.setInputTarget(cloud);
      align.setTargetFeatures(cloud_features);
      //Settings
      align.setMaximumIterations (5000); // Number of RANSAC iterations
      align.setNumberOfSamples (3); // Number of points to sample for generating/prerejecting a pose
      align.setCorrespondenceRandomness (3); // Number of nearest features to use
      align.setSimilarityThreshold (0.9f); // Polygonal edge length similarity threshold
      align.setMaxCorrespondenceDistance (2.5f * 0.01f);// Inlier threshold
      align.setInlierFraction (0.25f); // Required inlier fraction for accepting a pose hypothesis

      align.align(*cloud2_aligned);

      transformOut = align.getFinalTransformation ();

      if (align.hasConverged()){
        cout << " Successfully converged" << endl;
      }else {
        cout << "Alignment FAILED to converge" << endl;
      }
    }

    

  }



  
  // {
  //   pcl::ScopeTime t("Alignment");
  //   align.align (*cloud2_aligned);

  // }

  cout << "Completed alignment" << endl;

  // Eigen::Matrix4f transformOut;

  // transformOut = align.getFinalTransformation();

  cout << "Transform matrix out is: " << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    visu.addPointCloud (cloud, ColorHandlerT (cloud, 0.0, 255.0, 0.0), "cloud1");
    visu.addPointCloud (cloud2_aligned, ColorHandlerT (cloud2_aligned, 0.0, 0.0, 255.0), "cloud2_aligned");
    visu.addPointCloud (cloud2, ColorHandlerT (cloud2, 255.0, 0.0, 0.0), "cloud2");
    visu.spin ();
  }

  return transformOut;
}


void testLocalizationSequence(int argc, char** argv){
  // Input checks

  // Init
  int numberOfScans = 1;
  std::string filename;
  std::string filestem = "/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Data000";
  Eigen::Affine3f transform = Eigen::Affine3f::Identity();
  Eigen::Affine3f alignmentTransform = Eigen::Affine3f::Identity();

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

  // Load state
  Eigen::Array<float,6,10> state = readPathTextFile("/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Path.txt",10);
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

  

  // Setup to run sequence of scans 
  pcl::PCDReader reader;

  int option = 0;
  if (argc > 1){  
    numberOfScans = atoi(argv[1]);
    if (argc > 2){
      option = atoi(argv[2]);
    }
  }
  

  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));
  pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC(new pcl::PointCloud<pcl::PointNormal>(mtSurf, msSurf, pcl::PointNormal()));
  
  // Generate point cloud from object
  mp.pointCloudFromObject3D(0, msSurf, mtSurf, mapObjPC);

  Eigen::Vector3f rpy;

  // Run sequence of scans 
  for (int i = 0; i < numberOfScans; i++){
    cout << "Processing Scan " << i << endl;

    // Get filename:
    if (i < 9){
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
    // Reduce
    mp.meshFromScan(cloudReduced, cloud);
    cout << "Reduced PC" << endl;

    // Transform
    pcl::transformPointCloud(*cloudReduced, *cloudTransformed, transform);

    cout << "Transformed PC" << endl;

    // Run alignment
    
    alignmentTransform.matrix() = runAlignmentPreRejective(mapObjPC, cloudTransformed,option);

    // Update transform estimate - will need to check that this works
    transform = alignmentTransform * transform ;// alignmentTransform.transpose()

    cout << "State transform is " << transform.matrix() << endl;

    // // Transform and check scan?
    // if (i == 2 && firstTwo){
    //   i = i - 1;
    //   firstTwo = false;
    // }
    // transform.translation();
    state(0,i) = transform.matrix()(0,3);
    state(1,i) = transform.matrix()(1,3);
    state(2,i) = transform.matrix()(2,3);

    rpy = transform.rotation().eulerAngles(2, 1, 0); // Check ordering

    state(3,i) = rpy(2);
    state(4,i) = rpy(1);
    state(5,i) = rpy(0);


  }

  cout << "output state is:\n" << state << endl;
}

void testLocalizationSequenceWithOptions(int argc, char** argv){
  // Input checks

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
  int msSurf = 95; // These will be modified in the options below. 
  int mtSurf = 195;

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
      msSurf = 195;
      mtSurf = 395; 
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

  // Init mapping class with settings
  Mapping3D mp;
  mp.numRowsDesired = 125;
  mp.numColsDesired = 125;
  mp.maxNanAllowed = 20;
  mp.removeNanBuffer = 2;// Was 3
  // mp.nCtrlDefault[0] = 15;
  // mp.nCtrlDefault[1] = 15;

  // mp.newRowColBuffer = 20; // How many non new points in a row or column are permissible

  // New extension method
  mp.useNonRectData = true;

  int option = 2; // Option for localisation method (0 - PCL, 1 - RANSAC IA, 2 - Prerejective RANSAC)
  if (argc > 3){
    option = atoi(argv[3]);
  }

  // ------------------------------------------------------
  // ------------------------------------------------------
  // Load the object
  mp.addObjectFromFile((outFilestem + "final_obj.obj").c_str());
  // ------------------------------------------------------
  // ------------------------------------------------------

  // Load state
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

    // Setup to run sequence of scans 
  pcl::PCDReader reader;

  // Init clouds
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>); 
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));
  pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC(new pcl::PointCloud<pcl::PointNormal>(mtSurf, msSurf, pcl::PointNormal()));
  
  // Generate point cloud from object - to use for localisation
  mp.pointCloudFromObject3D(0, msSurf, mtSurf, mapObjPC);

  bool showAlignment = true;

  Eigen::Vector3f rpy; // init vector to store output state

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

    // Reduce
    mp.meshFromScan(cloudReduced, cloud);
    cout << "Reduced PC" << endl;

    // Transform
    pcl::transformPointCloud(*cloudReduced, *cloudTransformed, transform);
    cout << "Transformed PC" << endl;

    // Run alignment
    alignmentTransform.matrix() = runAlignmentPreRejective(mapObjPC, cloudTransformed,option,showAlignment);

    // Update transform estimate
    transform = alignmentTransform * transform ;

    cout << "State transform is " << transform.matrix() << endl;

    // Update State for tracking
    state(0,i) = transform.matrix()(0,3);
    state(1,i) = transform.matrix()(1,3);
    state(2,i) = transform.matrix()(2,3);

    rpy = transform.rotation().eulerAngles(2, 1, 0); // Check ordering

    state(3,i) = rpy(2);
    state(4,i) = rpy(1);
    state(5,i) = rpy(0);


  }

  // Basic for the moment - copy and paste state to get in Matlab - get an idea of tracking
  cout << "output state is:\n" << state << endl;

  // TDBM later - get last scan saved in PCD once aligned - how far off is it as a measure of localisation accuracy...
}

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "planar_seg_testing");
  ros::NodeHandle nh;

  // getAlignedCorrespondences(argc, argv);

  // testLocalizationSequence(argc, argv);

  testLocalizationSequenceWithOptions(argc, argv);

  // if (argc == 4){
  //   computePointCorrespondences(argc, argv);
    
  // }else{
  //   alignScans(argc, argv);
  // }
  
 
  // Spin
  // ros::spin ();
}
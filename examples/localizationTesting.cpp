#include <ros/ros.h>

//PCL includes
#include <Eigen/Core>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/time.h>

#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/fpfh_omp.h>

#include <pcl/common/transforms.h>

#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/sample_consensus_prerejective.h>
#include <pcl/registration/correspondence_estimation_normal_shooting.h>
#include <pcl/registration/correspondence_rejection_distance.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/visualization/pcl_visualizer.h>

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




int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "planar_seg_testing");
  ros::NodeHandle nh;

  if (argc == 4){
    computePointCorrespondences(argc, argv);
  }else{
    alignScans(argc, argv);
  }
  
 
  // Spin
  // ros::spin ();
}
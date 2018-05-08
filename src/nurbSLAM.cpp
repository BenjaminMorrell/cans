#include "cans/nurbSLAM.h"

using namespace PLib;

//-------------------------------------------------------------------
/*! 
  \brief  Default Constructor

  \warning 

  \author Benjamin Morrell
  \date 30 April 2018
*/
NurbSLAM::NurbSLAM(): 
    localisationOption(2), 
    alignmentOption(0),
    keypointOption(0),
    bRejectNonOverlappingInAlign(false),
    maxDistanceOverlap(0.2),
    bShowAlignment(false),
    bMapUpdatedFromScan(false), bMappingModeOn(false), bLocalisationModeOn(false),
    nSurfPointsFactor(5.0), 
    pclNormalRadiusSetting(0.05), pclFeatureRadiusSetting(0.1),
    ransac_inlierMultiplier(0.1), validInlierTheshold(0.5),inlierFraction(1.0), 
    modelResolutionKeypoints(0.005), minNeighboursKeypoints(5),
    ransac_maximumIterations(5000), ransac_numberOfSamples(3),
    ransac_correspondenceRandomness(3), ransac_similarityThreshold(0.9), ransac_inlierFraction(0.25),
    processTimes(5)
{
  state = Eigen::Affine3f::Identity();
  transformDelta = Eigen::Affine3f::Identity();

  // Mapping settings
  mp.numRowsDesired = 95;
  mp.numColsDesired = 95;
  mp.maxNanAllowed = 10;
  mp.removeNanBuffer = 2;
  mp.nCtrlDefault[0] = 17;
  mp.nCtrlDefault[1] = 17;
  mp.newRowColBuffer = 10; // How many non new points in a row or column are permissible
  // New extension method
  mp.useNonRectData = true;

  // nSurfPointsFactor = (float)mp.numRowsDesired/(float)mp.nCtrlDefault;

}

//-------------------------------------------------------------------
/*! 
  \brief  Destructor

  \author Benjamin Morrell
  \date 30 April 2018
*/
NurbSLAM::~NurbSLAM(){;}


//-------------------------------------------------------------------
// ------------- HGHER LEVEL FUNCTIONALITY  -----------------------------
//-------------------------------------------------------------------
void NurbSLAM::processScans(std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds){

  objectMeshList.clear();
  objIDList.clear();
  transformationList.clear();
  processTimes.clear();
  for (int i = 0; i < 5; i++){processTimes.push_back(0.0);}// reset to zero
  bMapUpdatedFromScan = false; // reset to false

  float msSurf;
  float mtSurf;

  std::chrono::high_resolution_clock::time_point t1;
  std::chrono::high_resolution_clock::time_point t2;
  std::chrono::duration<double, std::milli> dur;

  
  // Initial Scan processing
  for (int i = 0; i < clouds.size(); i++){

    objectMeshList.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired)));

    // Process the scan and perform data association
    objIDList.push_back(processSingleScan(clouds[i], objectMeshList[i]));

    
    cout << "\n\n\t\t FINISHED PROCESSING SCAN " << i << "\n\n";

    if (objIDList[i] == -9){
      cout << "Rejecting scan because of too many Nans (in NurbSLAM::processScans" << endl;
      continue;
    }

    if (!bMappingModeOn){
      t1 = std::chrono::high_resolution_clock::now();
      // Compute alignment for all matches
      if (objIDList[i] != -1){
        // Compute the loclisation transform for that object
        if (alignmentOption == 0){
          // Dense to dense alignment
          transformationList.push_back(alignScanWithMapObject(objIDList[i], objectMeshList[i]));
        }else if (alignmentOption == 1){
          // Keypoint to dense alignment
          transformationList.push_back(alignScanKeypointsWithMapObjectDense(objIDList[i], objectMeshList[i]));
        }else{
          // Keypoint to keypoint
          transformationList.push_back(alignScanKeypointsWithMapObjectKeypoints(objIDList[i], objectMeshList[i]));
        }

        if (inlierFraction < validInlierTheshold){
          cout << "inlier fraction of " << inlierFraction << " is below valid threshold: " << validInlierTheshold << ". Ignoring match." << endl;
          // Reject match - set to ignore flag
          objIDList[i] = -1;
          // Remove transformation
          transformationList.pop_back();
        }
      }

      t2 = std::chrono::high_resolution_clock::now();

      dur = t2-t1;
      processTimes[2] += dur.count();
      cout << "Duration for alignment is: " << dur.count() << "ms" << endl;

    }else{
      cout << "Mapping mode on, no alignment performed" << endl;
      processTimes[2] = -1.0; // Flag that there is no time
    }
  }

  cout << "Starting filter update" << endl;  
  // Perform the SLAM update with the computed transformations
  if (transformationList.size() > 0){
    t1 = std::chrono::high_resolution_clock::now();
    

    updateSLAMFilter(); // updates state


    t2 = std::chrono::high_resolution_clock::now();
    dur = t2-t1;
    processTimes[3] += dur.count();
    cout << "Duration for SLAM update filter is: " << dur.count() << "ms" << endl;

  }else{
    processTimes[3] = -1.0;
  }
  cout << "Updated filter" << endl;
  
  if (!bLocalisationModeOn){
    // Start time
    t1 = std::chrono::high_resolution_clock::now();
    

    // Update the map
    alignAndUpdateMeshes(); // uses the global lists


    cout << "Updated mesh" << endl;

    // Timing
    t2 = std::chrono::high_resolution_clock::now();
    dur = t2-t1;
    processTimes[4] += dur.count();
    cout << "Duration for Map update is: " << dur.count() << "ms" << endl;
  }else{
    processTimes[4] = -1.0;
  }
  
}


/*! 
  \brief Process scan for a single object

  mesh from scan
  adjust size
  transform
  search metrics
  data association

  \param[in]  cloud - input observation 
  \param[out] cloudTransformed - the output cloud - will be stored in an array objectMeshList

  \author Benjamin Morrell
  \date 2 May 2018
*/
int NurbSLAM::processSingleScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed){
  // cloudTransformed is what is the output

  // Initialise
  // clouds
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));
  // Search parameters
  std::vector<float> searchMetrics;
  int objID; 

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  // Process Scan to get mesh
  mp.meshFromScan(cloudReduced, cloud);

  if (mp.bRejectScan){
    cout << "Rejecting scan because of too many Nans" << endl;
    return -9;
  }

  // Resize cloudTransformed from size of cloudReduced
  if (cloudReduced->height < mp.numRowsDesired){
    pcl::common::deleteRows(*cloudTransformed, *cloudTransformed, std::max(1,(int)(1 + mp.numRowsDesired - cloudReduced->height)/2));
  }
  if (cloudReduced->width < mp.numColsDesired){
    pcl::common::deleteCols(*cloudTransformed, *cloudTransformed, std::max(1,(int)(1 + mp.numColsDesired - cloudReduced->width)/2));
  }

  // Transform with current state
  pcl::transformPointCloud(*cloudReduced, *cloudTransformed, state);

  cout << "Transformed point cloud" << endl;

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> dur = t2-t1;
  processTimes[0] += dur.count();
  cout << "Duration for mesh processing is: " << dur.count() << "ms" << endl;

  cout << "transform is: " << state.matrix() << endl;

  // Compute Metrics
  searchMetrics = mp.computeSearchMetrics(cloudTransformed);

  cout << "Computed Search Metrics:" << endl;
  for (int i = 0; i < 7; i++){cout << searchMetrics[i] << ", ";}cout << endl;

  // Data Association
  objID = mp.dataAssociation(searchMetrics);

  cout << "Completed data association" << endl;

  t1 = std::chrono::high_resolution_clock::now();
  dur = t1-t2;
  processTimes[1] += dur.count();
  cout << "Duration for data association is: " << dur.count() << "ms" << endl;

  return objID;
}


//-------------------------------------------------------------------
// ------------- LOCALISATION  -----------------------------
//-------------------------------------------------------------------
/*! 
  \brief  Compute the transform to match the scan to the map object using keypoints

  \author Benjamin Morrell
  \date 30 April 2018
*/
Eigen::Matrix4f NurbSLAM::alignScanKeypointsWithMapObjectKeypoints(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC){

  pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC_keypoints (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_keypoints (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr mapObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr obsObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_aligned (new pcl::PointCloud<pcl::PointNormal>);

  // Not use if this is needed
  // mapObjPC->is_dense = false;
  // obsObjPC->is_dense = false;

  std::cerr << "PointCloud1 dimensions, W: " << mapMeshList[objID]->width << "\tH: " << mapMeshList[objID]->height << "\t is dense? " << mapMeshList[objID]->is_dense << std::endl;
  std::cerr << "PointCloud2 dimensions, W: " << obsObjPC->width << "\tH: " << obsObjPC->height << "\t is dense? " << obsObjPC->is_dense << std::endl;

  cout << "Map obj at (0,0) is: " << mapMeshList[objID]->at(0,0) << endl;
  cout << "Obs obj at (0,0) is: " << obsObjPC->at(0,0) << endl;

  Eigen::Matrix4f transformOut;

  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  // Extract keypoints
  computeKeypoints(obsObjPC, search_method_, obsObjPC_keypoints);
  cout << "computed keypoints for cloud 1. have " << mapObjPC_keypoints->size() << endl;
  computeKeypoints(mapMeshList[objID], search_method_, mapObjPC_keypoints);
  cout << "computed keypoints for cloud 2. have " << obsObjPC_keypoints->size() << endl;
  // pcl::ISSKeypoint3D<pcl::PointNormal, pcl::PointNormal> iss_detector;

  // iss_detector.setSearchMethod(search_method_);
  // iss_detector.setSalientRadius (6 * modelResolutionKeypoints);
  // iss_detector.setNonMaxRadius (4 * modelResolutionKeypoints);
  // iss_detector.setThreshold21 (0.975);
  // iss_detector.setThreshold32 (0.975);
  // iss_detector.setMinNeighbors (minNeighboursKeypoints);
  // iss_detector.setNumberOfThreads (4);
  // iss_detector.setInputCloud (mapMeshList[objID]);
  // iss_detector.compute (*mapObjPC_keypoints);
  
  // iss_detector.setInputCloud (obsObjPC);
  // iss_detector.compute (*obsObjPC_keypoints);
  
  
  // Estimate normals
  pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(pclNormalRadiusSetting);
  // nest.setKSearch(7);
  nest.setSearchMethod(search_method_);
  nest.setInputCloud(obsObjPC);
  nest.compute(*obsObjPC);
  cout << "computed normals for cloud observation" << endl;
  
  // Estimate features
  pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  fest.setRadiusSearch(pclFeatureRadiusSetting);
  // fest.setKSearch(7);
  fest.setSearchMethod(search_method_);
  fest.setSearchSurface(mapMeshList[objID]);
  fest.setInputCloud(mapObjPC_keypoints);
  fest.setInputNormals(mapMeshList[objID]);
  fest.compute (*mapObjPC_features);
  fest.setSearchSurface(obsObjPC);
  fest.setInputCloud(obsObjPC_keypoints);
  fest.setInputNormals(obsObjPC);
  fest.compute (*obsObjPC_features);
  cout << "Features have been computed" << endl;
  

  //Perform alignment
  // Prerejective RANSAC
  // SampleConsensusPrerejective_Exposed<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  pcl::SampleConsensusPrerejective<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  align.setInputSource(obsObjPC_keypoints);
  align.setSourceFeatures(obsObjPC_features);
  align.setInputTarget(mapObjPC_keypoints);
  align.setTargetFeatures(mapObjPC_features);
  //Settings
  align.setMaximumIterations (ransac_maximumIterations); // Number of RANSAC iterations
  align.setNumberOfSamples (ransac_numberOfSamples); // Number of points to sample for generating/prerejecting a pose
  align.setCorrespondenceRandomness (ransac_correspondenceRandomness); // Number of nearest features to use
  align.setSimilarityThreshold (ransac_similarityThreshold); // Polygonal edge length similarity threshold
  align.setMaxCorrespondenceDistance (2.5f * ransac_inlierMultiplier);// Inlier threshold
  align.setInlierFraction (ransac_inlierFraction); // Required inlier fraction for accepting a pose hypothesis

  align.align(*obsObjPC_aligned);

  transformOut = align.getFinalTransformation ();

  if (align.hasConverged()){
    cout << "Successfully converged" << endl;
  }else {
    cout << "Alignment FAILED to converge" << endl;
  }
  

  cout << "Completed alignment" << endl;
  cout << "Transform matrix out is: " << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    visu.addPointCloud (mapMeshList[objID], ColorHandlerT (mapMeshList[objID], 0.0, 255.0, 0.0), "mapObjPC");
    visu.addPointCloud (obsObjPC_aligned, ColorHandlerT (obsObjPC_aligned, 0.0, 0.0, 255.0), "obsObjPC_alignedKP");
    visu.addPointCloud (obsObjPC, ColorHandlerT (obsObjPC, 255.0, 0.0, 0.0), "obsObjPC");
    visu.addPointCloud (mapObjPC_keypoints, ColorHandlerT (mapObjPC_keypoints, 0.0, 255.0, 255.0), "mapObjPCKP");
    visu.addPointCloud (obsObjPC_keypoints, ColorHandlerT (obsObjPC_keypoints, 255.0, 0.0, 255.0), "obsObjPCKP");
    visu.spin ();
  }

  return transformOut;
}

/*! 
  \brief  Compute the transform to match the scan to the map object using keypoints only from the obsection

  \author Benjamin Morrell
  \date 30 April 2018
*/
Eigen::Matrix4f NurbSLAM::alignScanKeypointsWithMapObjectDense(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC){

  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_keypoints (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr obsObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_aligned (new pcl::PointCloud<pcl::PointNormal>);

  // Not use if this is needed
  // mapObjPC->is_dense = false;
  // obsObjPC->is_dense = false;

  std::cerr << "PointCloud1 dimensions, W: " << mapMeshList[objID]->width << "\tH: " << mapMeshList[objID]->height << "\t is dense? " << mapMeshList[objID]->is_dense << std::endl;
  std::cerr << "PointCloud2 dimensions, W: " << obsObjPC->width << "\tH: " << obsObjPC->height << "\t is dense? " << obsObjPC->is_dense << std::endl;

  cout << "Map obj at (0,0) is: " << mapMeshList[objID]->at(0,0) << endl;
  cout << "Obs obj at (0,0) is: " << obsObjPC->at(0,0) << endl;

  Eigen::Matrix4f transformOut;

  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  // Extract keypoints
  computeKeypoints(obsObjPC, search_method_, obsObjPC_keypoints);

  // pcl::ISSKeypoint3D<pcl::PointNormal, pcl::PointNormal> iss_detector;

  // iss_detector.setSearchMethod(search_method_);
  // iss_detector.setSalientRadius (6 * modelResolutionKeypoints);
  // iss_detector.setNonMaxRadius (4 * modelResolutionKeypoints);
  // iss_detector.setThreshold21 (0.975);
  // iss_detector.setThreshold32 (0.975);
  // iss_detector.setMinNeighbors (minNeighboursKeypoints);
  // iss_detector.setNumberOfThreads (4);
  // iss_detector.setInputCloud (obsObjPC);
  // iss_detector.compute (*obsObjPC_keypoints);
  cout << "computed keypoints for observation. have " << obsObjPC_keypoints->size() << endl;
    
  // Estimate normals
  pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(pclNormalRadiusSetting);
  // nest.setKSearch(7);
  nest.setSearchMethod(search_method_);
  nest.setInputCloud(obsObjPC);
  nest.compute(*obsObjPC);
  cout << "computed normals for cloud observation" << endl;
  
  // Estimate features
  pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  fest.setRadiusSearch(pclFeatureRadiusSetting);
  // fest.setKSearch(7);
  fest.setSearchMethod(search_method_);
  fest.setSearchSurface(obsObjPC);
  fest.setInputCloud(obsObjPC_keypoints);
  fest.setInputNormals(obsObjPC);
  fest.compute (*obsObjPC_features);
  cout << "Features have been computed" << endl;
  

  //Perform alignment
  // Prerejective RANSAC
  // SampleConsensusPrerejective_Exposed<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  pcl::SampleConsensusPrerejective<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
  align.setInputSource(obsObjPC_keypoints);
  align.setSourceFeatures(obsObjPC_features);
  align.setInputTarget(mapMeshList[objID]);
  align.setTargetFeatures(mapFeatureList[objID]);
  //Settings
  align.setMaximumIterations (ransac_maximumIterations); // Number of RANSAC iterations
  align.setNumberOfSamples (ransac_numberOfSamples); // Number of points to sample for generating/prerejecting a pose
  align.setCorrespondenceRandomness (ransac_correspondenceRandomness); // Number of nearest features to use
  align.setSimilarityThreshold (ransac_similarityThreshold); // Polygonal edge length similarity threshold
  align.setMaxCorrespondenceDistance (2.5f * ransac_inlierMultiplier);// Inlier threshold
  align.setInlierFraction (ransac_inlierFraction); // Required inlier fraction for accepting a pose hypothesis

  align.align(*obsObjPC_aligned);

  transformOut = align.getFinalTransformation ();

  if (align.hasConverged()){
    cout << "Successfully converged" << endl;
  }else {
    cout << "Alignment FAILED to converge" << endl;
  }
  
  inlierFraction = (float)align.getInliers().size()/(float)(obsObjPC_keypoints->width*obsObjPC_keypoints->height);

  cout << "Completed alignment" << endl;
  cout << "Transform matrix out is: " << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    visu.addPointCloud (mapMeshList[objID], ColorHandlerT (mapMeshList[objID], 0.0, 255.0, 0.0), "mapObjPC");
    visu.addPointCloud (obsObjPC_aligned, ColorHandlerT (obsObjPC_aligned, 0.0, 0.0, 255.0), "obsObjPC_alignedKP");
    visu.addPointCloud (obsObjPC, ColorHandlerT (obsObjPC, 255.0, 0.0, 0.0), "obsObjPC");
    visu.addPointCloud (obsObjPC_keypoints, ColorHandlerT (obsObjPC_keypoints, 255.0, 0.0, 255.0), "obsObjPCKP");
    visu.spin ();
  }

  return transformOut;
}

/*! 
  \brief  Compute the transform to match the scan to the map object

  \author Benjamin Morrell
  \date 30 April 2018
*/
Eigen::Matrix4f NurbSLAM::alignScanWithMapObject(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC){
  
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr obsObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_aligned (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_filtered (new pcl::PointCloud<pcl::PointNormal>);

  // Not use if this is needed - is_dense means that there are no nans
  // mapObjPC->is_dense = false;
  // obsObjPC->is_dense = false;

  std::cerr << "PointCloud1 dimensions, W: " << mapMeshList[objID]->width << "\tH: " << mapMeshList[objID]->height << "\t is dense? " << mapMeshList[objID]->is_dense << std::endl;
  std::cerr << "PointCloud2 dimensions, W: " << obsObjPC->width << "\tH: " << obsObjPC->height << "\t is dense? " << obsObjPC->is_dense << std::endl;

  cout << "Map obj at (0,0) is: " << mapMeshList[objID]->at(0,0) << endl;
  cout << "Obs obj at (0,0) is: " << obsObjPC->at(0,0) << endl;

  Eigen::Matrix4f transformOut;

  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  if (localisationOption == 0){
    // ICP
    cout << "Starting alignment" << endl;
    pcl::IterativeClosestPoint<pcl::PointNormal, pcl::PointNormal> icp;
    // Set the input source and target
    icp.setInputSource (obsObjPC);
    icp.setInputTarget (mapMeshList[objID]);
    // Set the max correspondence distance to 5cm (e.g., correspondences with higher distances will be ignored)
    icp.setMaxCorrespondenceDistance (0.5);
    // Set the maximum number of iterations (criterion 1)
    icp.setMaximumIterations (ransac_maximumIterations);
    // Set the transformation epsilon (criterion 2)
    icp.setTransformationEpsilon (1e-8);
    // Set the euclidean distance difference epsilon (criterion 3)
    icp.setEuclideanFitnessEpsilon (1);
    
    // Perform the alignment
    icp.align (*obsObjPC_aligned);
    
    // Obtain the transformation that aligned cloud_source to cloud_source_registered
    transformOut = icp.getFinalTransformation ();
  }else{
    
    // Estimate normals
    pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
    nest.setRadiusSearch(pclNormalRadiusSetting);
    // nest.setKSearch(7);
    nest.setSearchMethod(search_method_);
    nest.setInputCloud(obsObjPC);
    nest.compute(*obsObjPC);
    cout << "computed normals for observation" << endl;

    cout << "Normals have been estimated" << endl;


    if (bRejectNonOverlappingInAlign){
      rejectNonOverlappingPoints( mapMeshList[objID], obsObjPC, obsObjPC_filtered);
    }else{
      obsObjPC_filtered = obsObjPC;
    }


    // Estimate features
    pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
    fest.setRadiusSearch(pclFeatureRadiusSetting);
    // fest.setKSearch(7);
    fest.setSearchMethod(search_method_);
    fest.setInputCloud(obsObjPC_filtered);
    fest.setSearchSurface(obsObjPC);
    fest.setInputNormals(obsObjPC);
    fest.compute (*obsObjPC_features);
    cout << "Features have been computed for observation" << endl;
  

    //Perform alignment 
    if (localisationOption == 1){
      // RANSAC Initial Alignment
      pcl::SampleConsensusInitialAlignment<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
      align.setInputSource(obsObjPC_filtered);
      align.setSourceFeatures(obsObjPC_features);
      align.setInputTarget(mapMeshList[objID]);
      align.setTargetFeatures(mapFeatureList[objID]);
      //Settings
      align.setMaximumIterations (ransac_maximumIterations); // Number of RANSAC iterations (1000 is default)
      align.setNumberOfSamples (ransac_numberOfSamples); // Number of points to sample for generating/prerejecting a pose (3 is default)
      align.setCorrespondenceRandomness (20); // Number of nearest features to use (default is 10)
      align.setMaxCorrespondenceDistance (0.1f);// Inlier threshold

      align.align(*obsObjPC_aligned);

      transformOut = align.getFinalTransformation ();

      if (align.hasConverged()){
        cout << " Successfully converged" << endl;
      }else {
        cout << "\n\tAlignment FAILED to converge\n" << endl;
      }

      

    }else if (localisationOption == 2){
      // Prerejective RANSAC
      // SampleConsensusPrerejective_Exposed<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
      pcl::SampleConsensusPrerejective<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
      align.setInputSource(obsObjPC_filtered);
      align.setSourceFeatures(obsObjPC_features);
      align.setInputTarget(mapMeshList[objID]);
      align.setTargetFeatures(mapFeatureList[objID]);
      //Settings
      align.setMaximumIterations (ransac_maximumIterations); // Number of RANSAC iterations
      align.setNumberOfSamples (ransac_numberOfSamples); // Number of points to sample for generating/prerejecting a pose
      align.setCorrespondenceRandomness (ransac_correspondenceRandomness); // Number of nearest features to use
      align.setSimilarityThreshold (ransac_similarityThreshold); // Polygonal edge length similarity threshold
      align.setMaxCorrespondenceDistance (2.5f * ransac_inlierMultiplier);// Inlier threshold
      align.setInlierFraction (ransac_inlierFraction); // Required inlier fraction for accepting a pose hypothesis

      align.align(*obsObjPC_aligned);

      transformOut = align.getFinalTransformation ();

      if (align.hasConverged()){
        cout << "Successfully converged" << endl;
      }else {
        cout << "Alignment FAILED to converge" << endl;
      }

      cout << "\n\n\t\t NUMBER OF INLIERS IS: " << align.getInliers().size() << "/ " << obsObjPC_filtered->width*obsObjPC_filtered->height << "\n\n\n";
      cout << "\n\n\t\t INLIER FRACTION IS: " << (float)align.getInliers().size()/(float)(obsObjPC_filtered->width*obsObjPC_filtered->height) << "\n\n\n";

      inlierFraction = (float)align.getInliers().size()/(float)(obsObjPC_filtered->width*obsObjPC_filtered->height);
    }
    
  }
    cout << "Completed alignment" << endl;
    cout << "Transform matrix out is: " << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    visu.addPointCloud (mapMeshList[objID], ColorHandlerT (mapMeshList[objID], 0.0, 255.0, 0.0), "mapObjPC");
    // visu.addPointCloudNormals (mapMeshList[objID], 100, 0.02, "mapObjPCNormals", 0);
    visu.addPointCloud (obsObjPC_aligned, ColorHandlerT (obsObjPC_aligned, 0.0, 0.0, 255.0), "obsObjPC_aligned");
    visu.addPointCloud (obsObjPC_filtered, ColorHandlerT (obsObjPC_filtered, 255.0, 0.0, 0.0), "obsObjPC");
    // visu.addPointCloudNormals (obsObjPC_filtered, 100, 0.02, "obsObjPCNormals", 0);
    visu.spin ();
  }

  return transformOut;
}

/*! 
  \brief  Compute the keypoints for a given scan

  \param cloud            input cloud
  \param search_method_   input search tree
  \param keypoints        output keypoints

  \author Benjamin Morrell
  \date 05 May 2018
*/
void NurbSLAM::computeKeypoints(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_, pcl::PointCloud<pcl::PointNormal>::Ptr keypoints){
   
  switch (keypointOption){
    case 0:
      {
      cout << "Computing Keypoints with an ISS detector" << endl;
      // Extract keypoints
      pcl::ISSKeypoint3D<pcl::PointNormal, pcl::PointNormal> iss_detector;
      iss_detector.setSearchMethod(search_method_);
      iss_detector.setSalientRadius (6 * modelResolutionKeypoints);
      iss_detector.setNonMaxRadius (4 * modelResolutionKeypoints);
      iss_detector.setThreshold21 (0.975);
      iss_detector.setThreshold32 (0.975);
      iss_detector.setMinNeighbors (minNeighboursKeypoints);
      iss_detector.setNumberOfThreads (4);
      iss_detector.setInputCloud (cloud);
      iss_detector.compute (*keypoints);
      }
      break;
    case 1:
      {
      cout << "NOT Computing Keypoints with a Harris3D detector - NOT IMPLEMENTED" << endl;
      // pcl::HarrisKeypoint6D<pcl::PointNormal, pcl::PointNormal, pcl::Normal> harris_detector;
      // harris_detector.setSearchMethod(search_method_);
      // harris_detector.setRadius(pclNormalRadiusSetting);
      // harris_detector.setNumberOfThreads (4);
      // // harris_detector.setNonMaxSupression(false);
      // // harris_detector.setThreshold(0.0);
      // harris_detector.setInputCloud (cloud);
      // // harris_detector.setNormals (cloud);
      // harris_detector.compute (*keypoints);
      }
      break;
    case 2:
      {
      cout << "Computing Keypoints with a Smoothed Surfaces detector" << endl;
      pcl::SmoothedSurfacesKeypoint<pcl::PointNormal, pcl::PointNormal> ss_detector;
      ss_detector.setSearchMethod(search_method_);
      ss_detector.setInputCloud(cloud);
      // ss_detector.setNormals(cloud);
      ss_detector.setNeighborhoodConstant(0.5);
      // ss_detector.setInputScale(0.0);
      }
      break;
    default :
      break;
  }

  // Others - harris6D...
  
  
}

/*! 
  \brief  Extract overlapping observation points 

  \author Benjamin Morrell
  \date 06 May 2018
*/
void NurbSLAM::rejectNonOverlappingPoints(pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsPCFilt){
  
  // EXPECT THERE TO BE NORMALS ALREADY ESTIMATED
  cout << "In rejectNonOverlappingPoints" << endl;
  pcl::CorrespondencesPtr correspondences (new pcl::Correspondences);
  pcl::CorrespondencesPtr corr_filtPtr (new pcl::Correspondences);

  // Estimate correspondences
  // bool doNormalShooting = true;

  // if (doNormalShooting){
  pcl::registration::CorrespondenceEstimationNormalShooting<pcl::PointNormal, pcl::PointNormal, pcl::PointNormal> corrEst;
  corrEst.setInputSource (obsObjPC);// Gets a corresponding point for every point in the Source
  corrEst.setSourceNormals (obsObjPC);
  corrEst.setInputTarget (mapObjPC); // Target for the source to find corresponding points in
  // Test the first 10 correspondences for each point in source, and return the best
  corrEst.setKSearch (10);
  cout << "Set up correspondences" << endl;
  corrEst.determineCorrespondences (*correspondences);
  // } else{
  //   pcl::registration::CorrespondenceEstimation<pcl::PointNormal, pcl::PointNormal, pcl::PointNormal> corrEst;
  //   corrEst.setInputSource (obsObjPC);// Gets a corresponding point for every point in the Source
  //   corrEst.setInputTarget (mapObjPC);
  //   corrEst.determineCorrespondences (*correspondences);
  // } 
  // #include <pcl/registration/correspondence_rejection_one_to_one.h>
  // pcl::registration::CorrespondenceRejectorOneToOne corrRej;
  // corrRej.applyRejection(*correspondences);

  cout << "estimated correspondences" << endl;

  // Correspondence rejection
  pcl::registration::CorrespondenceRejectorDistance corrRej;
  corrRej.setMaximumDistance( maxDistanceOverlap); 
  corrRej.setInputCorrespondences( correspondences);
  corrRej.getCorrespondences( *corr_filtPtr);

  pcl::Correspondences& corr_filt = *corr_filtPtr;

  // cout << "Corr filt is: " << corr_filt << endl;
  cout << "Rejected correspondences" << endl;

  // Get indices
  // std::vector<int> indices;
  pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());
  for (int i = 0; i < corr_filt.size(); i++){
    if (corr_filt[i].index_match >= 0){
      inliers->indices.push_back(corr_filt[i].index_query);       
    }
  }

  // cout << "Indices are: " << indices << endl;
  cout << "Made inliers object. Size of inliers is: " <<  inliers->indices.size() << endl;

  // Extract a cloud
  pcl::ExtractIndices<pcl::PointNormal> extract;
  extract.setInputCloud (obsObjPC);
  extract.setIndices (inliers);
  // extract.setIndices (corr_filt);
  extract.setNegative (false);
  extract.filter (*obsPCFilt);

  
  cout << "Extracted cloud. Size of output cloud is: " << obsPCFilt->width*obsPCFilt->height << endl;
}
//-------------------------------------------------------------------
// ------------- SLAM FUNCTIONS  -----------------------------
//-------------------------------------------------------------------
/*! 
  \brief  Perform the SLAM filter updates

  \warning currently a placeholder that just uses the first scan

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::updateSLAMFilter(){

  // PLACEHOLDER BEFORE A MORE IN-DEPTH SLAM ALGORITHM
  transformDelta.matrix() = transformationList[0];// list is of 4x4 matrices

  // Update the state with the computed transformation from the first observation
  state = transformDelta * state ;

  cout << "State transform is " << state.matrix() << endl;

}


/*! 
  \brief  Update the alignment of meshes and update the map

  the "update Map" step

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::alignAndUpdateMeshes(){
  
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));

  std::vector<float> searchMetrics;

  int updateID = -1;

  cout << "In Align and update Meshes. TransformDelta is:\n" << transformDelta.matrix() << endl;

  if (objectMeshList.size() > 0){
    cout << "Object mesh list has values cloud at [0] of size: " << objectMeshList[0]->size() << endl;
    cout << "and at (0,0) = " << objectMeshList[0]->at(0,0) << endl;
  }else{
    cout << "There are no objects to udpate" << endl;
  }

  for (int i = 0; i < objIDList.size(); i++){
    if (objIDList[i] == -2){
      // Ignore this scan
      cout << "ignoring scan due to bad alignment" << endl;
      continue;
    }else if (objIDList[i] == -9){
      // Ignore this scan
      cout << "ignoring scan due to too many Nans" << endl;
      continue;
    }
    // Align scans with new transformation delta (delat from existing transform)
    pcl::transformPointCloud(*objectMeshList[i], *cloudTransformed, transformDelta);

    if (objIDList[i] == -1){
      // New object
      cout << "\n\t\t ADDING NEW OBJECT \n\n";
      searchMetrics = mp.computeSearchMetrics(cloudTransformed);
      mp.addObject(cloudTransformed, searchMetrics);

      // ID for the latest object
      updateID = mp.objectMap.size() - 1;
    }else{
      cout << "\n\t\t UPDATING OBJECT " << objIDList[i] << "...\n\n";
      // Update object
      // Get size before update
      int nRows = mp.objectMap[objIDList[i]].ctrlPnts().rows();
      int nCols = mp.objectMap[objIDList[i]].ctrlPnts().cols();

      mp.updateObject(objIDList[i], cloudTransformed); 

      // track if there were any updates - change in the number of control points
      if (nRows != mp.objectMap[objIDList[i]].ctrlPnts().rows() ||
      nCols != mp.objectMap[objIDList[i]].ctrlPnts().cols()){
        cout << "Surface was updated" << endl;
        updateID = objIDList[i];
      }else{
        // Flag for no update
        updateID = -1;
      }
      
    }

    if (!bMappingModeOn){
      // Update mesh and feature for the items in the map
      updatePointCloudAndFeaturesInMap(updateID);
    }else{
      cout << "In Mapping mode. data points or map features for localisation to be computed" << endl;
      // Set flag that there have been updates
      bMapUpdatedFromScan = true;
    }
  }
}

/*! 
  \brief Computes the normals and features to store for the map object

  \param objID - object to update. -1 for no update

  Adds to or changes mapMeshlist and mapFeatureList

  \author Benjamin Morrell
  \date 2 May 2018
*/
void NurbSLAM::updatePointCloudAndFeaturesInMap(int objID){

  if (objID == -1){
    cout << "Nothing to update for normals and features" << endl;
    return;
  }

  // Get point cloud for object
  float msSurf = mp.objectMap[objID].ctrlPnts().rows()*nSurfPointsFactor;
  float mtSurf = mp.objectMap[objID].ctrlPnts().cols()*nSurfPointsFactor;
  cout << "surface points for object " << objID << " being updated with (ms, mt) = (" << msSurf << ", " << mtSurf << ")\n";

  // Initialise
  if (mapMeshList.size() <= objID){
    // New object
    cout << "new object - first mesh to generate" << endl;
    mapMeshList.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr (new pcl::PointCloud<pcl::PointNormal>(mtSurf, msSurf, pcl::PointNormal())));
    mapFeatureList.push_back(pcl::PointCloud<pcl::FPFHSignature33>::Ptr (new pcl::PointCloud<pcl::FPFHSignature33>));
  }else{
    // write a new object over the old
    cout << "update object - replace the mesh and features" << endl;
    mapMeshList[objID] = pcl::PointCloud<pcl::PointNormal>::Ptr (new pcl::PointCloud<pcl::PointNormal>(mtSurf, msSurf, pcl::PointNormal()));
    mapFeatureList[objID] = pcl::PointCloud<pcl::FPFHSignature33>::Ptr (new pcl::PointCloud<pcl::FPFHSignature33>);
  }
  
  // Generate point cloud for NURBS
  mp.pointCloudFromObject3D(objID, msSurf, mtSurf, mapMeshList[objID]);

  // Search tree
  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);

  // Estimate normals
  pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(pclNormalRadiusSetting);
  // nest.setKSearch(7);
  nest.setSearchMethod(search_method_);
  nest.setInputCloud(mapMeshList[objID]);
  nest.compute(*mapMeshList[objID]);
  cout << "computed normals for map cloud" << endl;

  // Estimate features
  pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  fest.setRadiusSearch(pclFeatureRadiusSetting);
  // fest.setKSearch(7);
  fest.setInputCloud(mapMeshList[objID]);
  fest.setSearchMethod(search_method_);
  fest.setInputNormals(mapMeshList[objID]);
  fest.compute (*mapFeatureList[objID]);
  cout << "Features have been computed for map cloud" << endl;

  // Set flag that there have been updates
  bMapUpdatedFromScan = true;

}

/*! 
  \brief  Sets the state
  
  Used for initialisation (if not starting at zero)
  Also used in mapping paradigms

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::setState(Eigen::Affine3f inputState){
  this->state = inputState;
}


/*! 
  \brief  Get the current state

  \author Benjamin Morrell
  \date 30 April 2018
*/
Eigen::Affine3f NurbSLAM::getState(){
  return state;
}

/*! 
  \brief Turns on mapping mode

  \warning only one mode is on at a time

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::activateMappingMode(){
  bMappingModeOn = true;
  bLocalisationModeOn = false;
}

/*! 
  \brief Turns on mapping mode

  \warning only one mode is on at a time
  \warning Need to load an object for this to work

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::activateLocalisationMode(){
  bLocalisationModeOn = true;
  bMappingModeOn = false;
}

/*! 
  \brief Loads object into the map

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::loadObjectIntoMap(std::string filename){
  
  // Add object to the map
  mp.addObjectFromFile(filename.c_str());

  // Get the object ID for the new object
  int objID = mp.objectMap.size() - 1;

  // Process the point clouds for the new object
  updatePointCloudAndFeaturesInMap(objID);

}


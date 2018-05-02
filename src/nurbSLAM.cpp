#include "cans/nurbSLAM.h"

using namespace PLib;

//-------------------------------------------------------------------
/*! 
  \brief  Default Constructor

  \warning 

  \author Benjamin Morrell
  \date 30 April 2018
*/
NurbSLAM::NurbSLAM(): localisationOption(2), bShowAlignment(false),
    nSurfPointsFactor(3.0), pclNormalRadiusSetting(0.05), pclFeatureRadiusSetting(0.1),
    bUseKeypoints(true), modelResolution(0.005), inlierMultiplier(0.1)
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

  float msSurf;
  float mtSurf;

  // Initial Scan processing
  for (int i = 0; i < clouds.size(); i++){

    objectMeshList.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired)));

    // Process the scan and perform data association
    objIDList.push_back(processSingleScan(clouds[i], objectMeshList[i]));

    cout << "\n\n\t\t FINISHED PROCESSING SCAN " << i << "\n\n";

    // Compute alignment for all matches
    if (objIDList[i] != -1){
      
      // Get point cloud for object
      msSurf = mp.objectMap[objIDList[i]].ctrlPnts().rows()*nSurfPointsFactor;
      mtSurf = mp.objectMap[objIDList[i]].ctrlPnts().cols()*nSurfPointsFactor;
      cout << "surface points for object " << objIDList[i] << " being updated with (ms, mt) = (" << msSurf << ", " << mtSurf << ")\n";
      pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC(new pcl::PointCloud<pcl::PointNormal>(mtSurf, msSurf, pcl::PointNormal()));
      mp.pointCloudFromObject3D(objIDList[i], msSurf, mtSurf, mapObjPC);

      cout << "Map obj at (0,0) is: " << mapObjPC->at(0,0) << endl;
      
      // Compute the loclisation transform for that object
      if (bUseKeypoints){
        transformationList.push_back(alignScanKeypointsWithMapObject(mapObjPC, objectMeshList[i]));
      }else{
        transformationList.push_back(alignScanWithMapObject(mapObjPC, objectMeshList[i]));
      }
    }
  }

  cout << "Starting filter update" << endl;
  // Perform the SLAM update with the computed transformations
  if (transformationList.size() > 0){
    updateSLAMFilter(); // updates state
  }

  cout << "Updated filter" << endl;

  // Update the map
  alignAndUpdateMeshes(); // uses the global lists

  cout << "Updated mesh" << endl;
  
}


// Currently for a single scan - multiple scans to come later
int NurbSLAM::processSingleScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed){
  
  // Initialise
  // clouds
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));
  // Search parameters
  std::vector<float> searchMetrics;
  int objID; 

  // Process Scan to get mesh
  mp.meshFromScan(cloudReduced, cloud);

  cout << "transform is: " << state.matrix() << endl;

  // Transform with current state
  pcl::transformPointCloud(*cloudReduced, *cloudTransformed, state);

  cout << "Transformed point cloud" << endl;

  // Compute Metrics
  searchMetrics = mp.computeSearchMetrics(cloudTransformed);

  cout << "Computed Search Metrics:" << endl;
  for (int i = 0; i < 7; i++){cout << searchMetrics[i] << ", ";}cout << endl;

  // Data Association
  objID = mp.dataAssociation(searchMetrics);

  cout << "Completed data association" << endl;

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
Eigen::Matrix4f NurbSLAM::alignScanKeypointsWithMapObject(pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC){
  pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC_keypoints (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_keypoints (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr mapObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr obsObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_aligned (new pcl::PointCloud<pcl::PointNormal>);

  // Not use if this is needed
  // mapObjPC->is_dense = false;
  // obsObjPC->is_dense = false;

  std::cerr << "PointCloud1 dimensions, W: " << mapObjPC->width << "\tH: " << mapObjPC->height << "\t is dense? " << mapObjPC->is_dense << std::endl;
  std::cerr << "PointCloud2 dimensions, W: " << obsObjPC->width << "\tH: " << obsObjPC->height << "\t is dense? " << obsObjPC->is_dense << std::endl;

  cout << "Map obj at (0,0) is: " << mapObjPC->at(0,0) << endl;
  cout << "Obs obj at (0,0) is: " << obsObjPC->at(0,0) << endl;

  Eigen::Matrix4f transformOut;

  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  // Extract keypoints
  pcl::ISSKeypoint3D<pcl::PointNormal, pcl::PointNormal> iss_detector;

  iss_detector.setSearchMethod(search_method_);
  iss_detector.setSalientRadius (6 * modelResolution);
  iss_detector.setNonMaxRadius (4 * modelResolution);
  iss_detector.setThreshold21 (0.975);
  iss_detector.setThreshold32 (0.975);
  iss_detector.setMinNeighbors (5);
  iss_detector.setNumberOfThreads (4);
  iss_detector.setInputCloud (mapObjPC);
  iss_detector.compute (*mapObjPC_keypoints);
  cout << "computed keypoints for cloud 1. have " << mapObjPC_keypoints->size() << endl;
  iss_detector.setInputCloud (obsObjPC);
  iss_detector.compute (*obsObjPC_keypoints);
  cout << "computed keypoints for cloud 2. have " << obsObjPC_keypoints->size() << endl;
    
  // Estimate normals
  pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(pclNormalRadiusSetting);
  // nest.setKSearch(7);
  nest.setSearchMethod(search_method_);
  nest.setInputCloud(mapObjPC);
  nest.compute(*mapObjPC);
  cout << "computed normals for cloud 1" << endl;
  nest.setInputCloud(obsObjPC);
  nest.compute(*obsObjPC);
  cout << "computed normals for cloud 2" << endl;

  cout << "Normals have been estimated" << endl;

  // Estimate features
  pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  fest.setRadiusSearch(pclFeatureRadiusSetting);
  // fest.setKSearch(7);
  fest.setSearchMethod(search_method_);
  fest.setSearchSurface(mapObjPC);
  fest.setInputCloud(mapObjPC_keypoints);
  fest.setInputNormals(mapObjPC);
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
  align.setMaximumIterations (5000); // Number of RANSAC iterations
  align.setNumberOfSamples (3); // Number of points to sample for generating/prerejecting a pose
  align.setCorrespondenceRandomness (3); // Number of nearest features to use
  align.setSimilarityThreshold (0.9f); // Polygonal edge length similarity threshold
  align.setMaxCorrespondenceDistance (2.5f * inlierMultiplier);// Inlier threshold
  align.setInlierFraction (0.25f); // Required inlier fraction for accepting a pose hypothesis

  align.align(*obsObjPC_aligned);

  transformOut = align.getFinalTransformation ();

  if (align.hasConverged()){
    cout << " Successfully converged" << endl;
  }else {
    cout << "Alignment FAILED to converge" << endl;
  }
  

  cout << "Completed alignment" << endl;
  cout << "Transform matrix out is: " << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    visu.addPointCloud (mapObjPC, ColorHandlerT (mapObjPC, 0.0, 255.0, 0.0), "mapObjPC");
    visu.addPointCloud (obsObjPC_aligned, ColorHandlerT (obsObjPC_aligned, 0.0, 0.0, 255.0), "obsObjPC_alignedKP");
    visu.addPointCloud (obsObjPC, ColorHandlerT (obsObjPC, 255.0, 0.0, 0.0), "obsObjPC");
    visu.addPointCloud (mapObjPC_keypoints, ColorHandlerT (mapObjPC_keypoints, 0.0, 255.0, 255.0), "mapObjPCKP");
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
Eigen::Matrix4f NurbSLAM::alignScanWithMapObject(pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC){
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr mapObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr obsObjPC_features (new pcl::PointCloud<pcl::FPFHSignature33>);
  pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC_aligned (new pcl::PointCloud<pcl::PointNormal>);

  // Not use if this is needed - is_dense means that there are no nans
  // mapObjPC->is_dense = false;
  // obsObjPC->is_dense = false;

  std::cerr << "PointCloud1 dimensions, W: " << mapObjPC->width << "\tH: " << mapObjPC->height << "\t is dense? " << mapObjPC->is_dense << std::endl;
  std::cerr << "PointCloud2 dimensions, W: " << obsObjPC->width << "\tH: " << obsObjPC->height << "\t is dense? " << obsObjPC->is_dense << std::endl;

  cout << "Map obj at (0,0) is: " << mapObjPC->at(0,0) << endl;
  cout << "Obs obj at (0,0) is: " << obsObjPC->at(0,0) << endl;

  Eigen::Matrix4f transformOut;

  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  if (localisationOption == 0){
    // ICP
    cout << "Starting alignment" << endl;
    pcl::IterativeClosestPoint<pcl::PointNormal, pcl::PointNormal> icp;
    // Set the input source and target
    icp.setInputSource (obsObjPC);
    icp.setInputTarget (mapObjPC);
    // Set the max correspondence distance to 5cm (e.g., correspondences with higher distances will be ignored)
    icp.setMaxCorrespondenceDistance (0.5);
    // Set the maximum number of iterations (criterion 1)
    icp.setMaximumIterations (500);
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
    nest.setInputCloud(mapObjPC);
    nest.compute(*mapObjPC);
    cout << "computed normals for cloud 1" << endl;
    nest.setInputCloud(obsObjPC);
    nest.compute(*obsObjPC);
    cout << "computed normals for cloud 2" << endl;

    cout << "Normals have been estimated" << endl;

    // Estimate features
    pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
    fest.setRadiusSearch(pclFeatureRadiusSetting);
    // fest.setKSearch(7);
    fest.setInputCloud(mapObjPC);
    fest.setSearchMethod(search_method_);
    fest.setInputNormals(mapObjPC);
    fest.compute (*mapObjPC_features);
    fest.setInputCloud(obsObjPC);
    fest.setInputNormals(obsObjPC);
    fest.compute (*obsObjPC_features);
    cout << "Features have been computed" << endl;
  

    //Perform alignment 
    if (localisationOption == 1){
      // RANSAC Initial Alignment
      pcl::SampleConsensusInitialAlignment<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> align;
      align.setInputSource(obsObjPC);
      align.setSourceFeatures(obsObjPC_features);
      align.setInputTarget(mapObjPC);
      align.setTargetFeatures(mapObjPC_features);
      //Settings
      align.setMaximumIterations (500); // Number of RANSAC iterations (1000 is default)
      align.setNumberOfSamples (3); // Number of points to sample for generating/prerejecting a pose (3 is default)
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
      align.setInputSource(obsObjPC);
      align.setSourceFeatures(obsObjPC_features);
      align.setInputTarget(mapObjPC);
      align.setTargetFeatures(mapObjPC_features);
      //Settings
      align.setMaximumIterations (5000); // Number of RANSAC iterations
      align.setNumberOfSamples (3); // Number of points to sample for generating/prerejecting a pose
      align.setCorrespondenceRandomness (3); // Number of nearest features to use
      align.setSimilarityThreshold (0.9f); // Polygonal edge length similarity threshold
      align.setMaxCorrespondenceDistance (2.5f * 0.01f);// Inlier threshold
      align.setInlierFraction (0.25f); // Required inlier fraction for accepting a pose hypothesis

      align.align(*obsObjPC_aligned);

      transformOut = align.getFinalTransformation ();

      if (align.hasConverged()){
        cout << " Successfully converged" << endl;
      }else {
        cout << "Alignment FAILED to converge" << endl;
      }
    }

    cout << "Completed alignment" << endl;
    cout << "Transform matrix out is: " << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    visu.addPointCloud (mapObjPC, ColorHandlerT (mapObjPC, 0.0, 255.0, 0.0), "mapObjPC");
    visu.addPointCloud (obsObjPC_aligned, ColorHandlerT (obsObjPC_aligned, 0.0, 0.0, 255.0), "obsObjPC_aligned");
    visu.addPointCloud (obsObjPC, ColorHandlerT (obsObjPC, 255.0, 0.0, 0.0), "obsObjPC");
    visu.spin ();
  }

  return transformOut;

  }
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

  cout << "In Align and update Meshes. TransformDelta is: " << transformDelta.matrix() << endl;

  cout << "Object mesh list has values cloud at [0] of size: " << objectMeshList[0]->size() << endl;
  cout << "and at (0,0) = " << objectMeshList[0]->at(0,0) << endl;

  for (int i = 0; i < objIDList.size(); i++){
    // Align scans with new transformation delta (delat from existing transform)
    pcl::transformPointCloud(*objectMeshList[i], *cloudTransformed, transformDelta);

    if (objIDList[i] == -1){
      // New object
      cout << "\n\t\t ADDING NEW OBJECT \n\n";
      searchMetrics = mp.computeSearchMetrics(cloudTransformed);
      mp.addObject(cloudTransformed, searchMetrics);
    }else{
      cout << "\n\t\t UPDATING OBJECT " << objIDList[i] << "...\n\n";
      // Update object
      mp.updateObject(objIDList[i], cloudTransformed); 
    }
  }

}

/*! 
  \brief  Initialise the state (if not starting at zero)

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::initState(Eigen::Affine3f startingState){
  this->state = startingState;
}

/*! 
  \brief  Get the current state

  \author Benjamin Morrell
  \date 30 April 2018
*/
Eigen::Affine3f NurbSLAM::getState(){
  return state;
}

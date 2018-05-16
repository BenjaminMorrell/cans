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
    processTimes(5), bObjectNormalsComputed(false), processModel(0),numberOfPointsInAlignment(0),
    bUseFullAlignmentTransformInUpdate(false), bRejectAlignment(false), bUseOldStateForNewObjects(false),
    rejectCriteria(6), bKeepPConstant(false), mapCountThreshold(3), mapExtendThreshold(2), updateCount(0),
    bLocalisationRejectionOn(false),lastMatchID(-1)
{// TODO - make this initialiser list organised...
  state = Eigen::Affine3f::Identity();
  oldState = Eigen::Affine3f::Identity();
  previousState = Eigen::Affine3f::Identity();
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
  cout << "Setting up EKF matrices" << endl;
  // EKF init
  pNoisePos = 0.05;
  pNoiseVel = 0.05;
  pNoiseAccel = 0.03;
  pNoiseAng = 0.01;
  pNoiseMultiplier = 1.0;
  qNoiseMultiplier = 3.0;

  // Observation noise
  noiseObsBasePos = 1.0;
  noiseObsMultPos = 2.0;
  noiseObsMultPosErr = 1.0;
  noiseObsBaseAng = 0.5;
  noiseObsMultAng = 0.5;
  noiseObsMultAngErr = 0.25;

  rMatMultiplier = 1.0;

  // Set up matrices
  setInitEKFStates();

  rejectCriteria[0] = 0.78; // 45 degrees
  rejectCriteria[1] = 1.57;
  rejectCriteria[2] = 1.5; // Linear error
  rejectCriteria[3] = 3.0;
  rejectCriteria[4] = 0.6; // inlier criteria
  rejectCriteria[5] = 0.25; // fraction of points compared to desired. 

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
void NurbSLAM::processScans(std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds, float timestep){

  objectMeshList.clear();
  objIDList.clear();
  transformationList.clear();
  inlierFractionList.clear();
  processTimes.clear();
  transformDelta = Eigen::Affine3f::Identity();

  // Store the previous state
  previousState = state;

  for (int i = 0; i < 5; i++){processTimes.push_back(0.0);}// reset to zero
  bMapUpdatedFromScan = false; // reset to false

  float msSurf;
  float mtSurf;

  std::chrono::high_resolution_clock::time_point t1;
  std::chrono::high_resolution_clock::time_point t2;
  std::chrono::duration<double, std::milli> dur;

  // Process step of EKF
  processStepEKF(timestep);

  oldState = state;

  int preAlignmentID;
  
  // Initial Scan processing
  for (int i = 0; i < clouds.size(); i++){

    bObjectNormalsComputed = false;

    objectMeshList.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired)));

    // Process the scan and perform data association
    objIDList.push_back(processSingleScan(clouds[i], objectMeshList[i]));

    preAlignmentID = objIDList[i];
    
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

        inlierFractionList.push_back(inlierFraction);

        if (inlierFraction < validInlierTheshold){
          cout << "inlier fraction of " << inlierFraction << " is below valid threshold: " << validInlierTheshold << ". Ignoring match." << endl;
          // Reject match - set to ignore flag
          objIDList[i] = -1;
          // Remove transformation
          transformationList.pop_back();
          inlierFractionList.pop_back();

        }
      }

      if (bLocalisationModeOn && objIDList[i] == -1){
        // If doing localisation and there is a failure - try use the last matched object (if different)
        cout << "Failed alignment on obj: " << objIDList[i] << ", last match is: " << lastMatchID << endl;
        if (lastMatchID != preAlignmentID || lastMatchID != -1){
          // Try localisation on lastMatchID
          cout << "Trying alignment with previous match, of ID: " << lastMatchID << endl;
          transformationList.push_back(alignScanWithMapObject(lastMatchID, objectMeshList[i]));

          inlierFractionList.push_back(inlierFraction);

          if (inlierFraction < validInlierTheshold){
            cout << "inlier fraction of " << inlierFraction << " is below valid threshold: " << validInlierTheshold << ". Ignoring match." << endl;
            objIDList[i] = -1;
            // Remove transformation
            transformationList.pop_back();
            inlierFractionList.pop_back();
          }else{
            cout << "Success in alignment with previous match!" << endl;
            objIDList[i] = lastMatchID;
          }
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

  if (objIDList[0] != -1){
    lastMatchID = objIDList[0];
  }
  

  // Update map match count
  bool bIsMatched;
  cout << "Match count is: [";
  for (int i = 0; i < mapMatchCount.size(); i++){
    bIsMatched = false;
    for (int j = 0; j < objIDList.size(); j++){
      if (i == objIDList[j]){
        // Object has been matched, reset count to zero
        mapMatchCount[i] = 0;
        bIsMatched = true;
      }
    }
    if (!bIsMatched){
      // Increment if not matched
      mapMatchCount[i]++;
    }
    cout << mapMatchCount[i] << ", ";
  }
  cout << "]" << endl;
  

  cout << "Starting filter update" << endl;  
  // Perform the SLAM update with the computed transformations
  if (transformationList.size() > 0){
    t1 = std::chrono::high_resolution_clock::now();
    

    updateSLAMFilter(timestep); // updates state and EKF states


    t2 = std::chrono::high_resolution_clock::now();
    dur = t2-t1;
    processTimes[3] += dur.count();
    cout << "Duration for SLAM update filter is: " << dur.count() << "ms" << endl;

  }else{
    processTimes[3] = -1.0;
  }
  cout << "Updated filter" << endl;

  if (bLocalisationRejectionOn && bLocalisationModeOn && (bRejectAlignment || objIDList[0] < 0)){
    // If doing localisation and a scan is rejected - reset to the previous state
    cout << "In Localisation Mode. Not changing state this iteration..." << endl;
    state = previousState;
  }
  
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
    // Null reading for localisation
    processTimes[4] = -1.0;

    // reset
    bRejectAlignment = false;
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

  cout << "transform is: \n" << state.matrix() << endl;

  // Compute Metrics
  searchMetrics = mp.computeSearchMetrics(cloudTransformed);

  cout << "Computed Search Metrics:" << endl;
  for (int i = 0; i < 7; i++){cout << searchMetrics[i] << ", ";}cout << endl;

  // Data Association
  objID = mp.dataAssociation(searchMetrics);

  cout << "Completed data association" << endl;

  if (objID >= 0){
    if (mapExtendCount[objID] > mapExtendThreshold && bMappingModeOn){
      cout << "No extension for too many observations. Creating new object" << endl;
      objID = -1;
    }
  }

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
  computeNormals(obsObjPC);
  // pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
  // nest.setRadiusSearch(pclNormalRadiusSetting);
  // // nest.setKSearch(7);
  // nest.setSearchMethod(search_method_);
  // nest.setViewPoint(state.matrix()(0,3),state.matrix()(1,3),state.matrix()(2,3)); // Set the viewpoint to the current pose
  // nest.setInputCloud(obsObjPC);
  // nest.compute(*obsObjPC);
  cout << "computed normals for cloud observation" << endl;
  bObjectNormalsComputed = true;
  
  // Estimate features
  computeFeatures(mapObjPC_keypoints, mapMeshList[objID], mapObjPC_features);
  computeFeatures(obsObjPC_keypoints, obsObjPC, obsObjPC_features);
  // pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  // fest.setRadiusSearch(pclFeatureRadiusSetting);
  // // fest.setKSearch(7);
  // fest.setSearchMethod(search_method_);
  // fest.setSearchSurface(mapMeshList[objID]);
  // fest.setInputCloud(mapObjPC_keypoints);
  // fest.setInputNormals(mapMeshList[objID]);
  // fest.compute (*mapObjPC_features);
  // fest.setSearchSurface(obsObjPC);
  // fest.setInputCloud(obsObjPC_keypoints);
  // fest.setInputNormals(obsObjPC);
  // fest.compute (*obsObjPC_features);
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
  cout << "Transform matrix out is: \n" << transformOut << endl;
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
  computeNormals(obsObjPC);
  // pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
  // nest.setRadiusSearch(pclNormalRadiusSetting);
  // // nest.setKSearch(7);
  // nest.setSearchMethod(search_method_);
  // nest.setViewPoint(state.matrix()(0,3),state.matrix()(1,3),state.matrix()(2,3)); // Set the viewpoint to the current pose
  // nest.setInputCloud(obsObjPC);
  // nest.compute(*obsObjPC);
  cout << "computed normals for cloud observation" << endl;
  bObjectNormalsComputed = true;
  
  // Estimate features
  computeFeatures(obsObjPC_keypoints, obsObjPC, obsObjPC_features);
  // pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  // fest.setRadiusSearch(pclFeatureRadiusSetting);
  // // fest.setKSearch(7);
  // fest.setSearchMethod(search_method_);
  // fest.setSearchSurface(obsObjPC);
  // fest.setInputCloud(obsObjPC_keypoints);
  // fest.setInputNormals(obsObjPC);
  // fest.compute (*obsObjPC_features);
  cout << "Features have been computed for observation keypoints" << endl;
  

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
  cout << "Transform matrix out is: \n" << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    // visu.setBackgroundColor (0, 0, 0);
    visu.addPointCloud (mapMeshList[objID], ColorHandlerT (mapMeshList[objID], 0.0, 255.0, 0.0), "mapObjPC");
    visu.addPointCloudNormals<pcl::PointNormal> (mapMeshList[objID], 50, 1.0, "mapObjPCNormals", 0);
    visu.addPointCloud (obsObjPC_aligned, ColorHandlerT (obsObjPC_aligned, 0.0, 0.0, 255.0), "obsObjPC_alignedKP");
    visu.addPointCloud (obsObjPC, ColorHandlerT (obsObjPC, 255.0, 0.0, 0.0), "obsObjPC");
    visu.addPointCloud (obsObjPC_keypoints, ColorHandlerT (obsObjPC_keypoints, 255.0, 0.0, 255.0), "obsObjPCKP");
    visu.addPointCloudNormals<pcl::PointNormal> (obsObjPC_aligned, 10, 0.2, "obsObjPCNormals", 0);
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

  // pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
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
    // pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
    // nest.setRadiusSearch(pclNormalRadiusSetting);
    // // nest.setKSearch(7);
    // nest.setSearchMethod(search_method_);
    // nest.setViewPoint(state.matrix()(0,3),state.matrix()(1,3),state.matrix()(2,3)); // Set the viewpoint to the current pose
    // nest.setInputCloud(obsObjPC);
    // nest.compute(*obsObjPC);
    computeNormals(obsObjPC);
    cout << "computed normals for observation" << endl;

    if (mapMatchCount[objID] > mapCountThreshold || bLocalisationModeOn){
      computeNormals(mapMeshList[objID]);
      cout << "\nComputed normals for Map because mapCount was: " << mapMatchCount[objID] << endl;
    }


    bObjectNormalsComputed = true;

    if (bRejectNonOverlappingInAlign){
      rejectNonOverlappingPoints( mapMeshList[objID], obsObjPC, obsObjPC_filtered);

      if (obsObjPC_filtered->width*obsObjPC_filtered->height < 10){
        cout << "Not enough overlap points" << endl;
        inlierFraction = 0.0;
        return transformOut;
      }
    }else{
      obsObjPC_filtered = obsObjPC;
    }

    numberOfPointsInAlignment = obsObjPC_filtered->width*obsObjPC_filtered->height;

    // Estimate features
    computeFeatures(obsObjPC_filtered, obsObjPC, obsObjPC_features);
    // pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
    // fest.setRadiusSearch(pclFeatureRadiusSetting);
    // // fest.setKSearch(7);
    // fest.setSearchMethod(search_method_);
    // fest.setInputCloud(obsObjPC_filtered);
    // fest.setSearchSurface(obsObjPC);
    // fest.setInputNormals(obsObjPC);
    // fest.compute (*obsObjPC_features);
    cout << "Features have been computed for observation" << endl;

     if (mapMatchCount[objID] > mapCountThreshold || bLocalisationModeOn){
      computeFeatures(mapMeshList[objID], mapMeshList[objID], mapFeatureList[objID]);
      cout << "\nComputed features for Map because mapCount was: " << mapMatchCount[objID] << endl;
    }
  

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
    cout << "Transform matrix out is: \n" << transformOut << endl;
  if (bShowAlignment){
    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment");
    // visu.setBackgroundColor (0, 0, 0);
    visu.addPointCloud (mapMeshList[objID], ColorHandlerT (mapMeshList[objID], 0.0, 255.0, 0.0), "mapObjPC");
    visu.addPointCloudNormals<pcl::PointNormal> (mapMeshList[objID], 50, 1.0, "mapObjPCNormals", 0);
    visu.addPointCloud (obsObjPC_aligned, ColorHandlerT (obsObjPC_aligned, 0.0, 0.0, 255.0), "obsObjPC_aligned");
    visu.addPointCloud (obsObjPC_filtered, ColorHandlerT (obsObjPC_filtered, 255.0, 0.0, 0.0), "obsObjPC");
    visu.addPointCloudNormals<pcl::PointNormal, pcl::PointNormal> (obsObjPC_aligned, obsObjPC_aligned, 10, 0.2, "obsObjPCNormals", 0);
    visu.spin ();
  }

  // bool bSaveAlignmentData = false;

  // if (bSaveAlignmentData){

  // }

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
  \brief  Compute the normals for a given scan

  \param cloud            input cloud, and the cloud where the normals are stored 

  \author Benjamin Morrell
  \date 08 May 2018
*/
void NurbSLAM::computeNormals(pcl::PointCloud<pcl::PointNormal>::Ptr cloud){
  
  // Init search tree
  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  // Estimate normals
  pcl::NormalEstimationOMP<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(pclNormalRadiusSetting);
  // nest.setKSearch(7);
  nest.setSearchMethod(search_method_);
  nest.setViewPoint(state.matrix()(0,3),state.matrix()(1,3),state.matrix()(2,3)); // Set the viewpoint to the current pose
  nest.setInputCloud(cloud);
  nest.compute(*cloud);
  cout << "computed normals for cloud" << endl;
  
}

/*! 
  \brief  Compute the features for a given scan

  \param cloud            downsampled cloud
  \param searchSurface    original cloud for searching neighbours
  \param[out] features    the output features

  \warning It is expected that the normals are already computed

  \author Benjamin Morrell
  \date 08 May 2018
*/
void NurbSLAM::computeFeatures(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr searchSurface, pcl::PointCloud<pcl::FPFHSignature33>::Ptr features){
  
  // Init search tree
  pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);
  
  // Estimate features
  pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  fest.setRadiusSearch(pclFeatureRadiusSetting);
  // fest.setKSearch(7);
  fest.setSearchMethod(search_method_);
  fest.setInputCloud(cloud);
  fest.setSearchSurface(searchSurface);
  fest.setInputNormals(searchSurface);
  fest.compute (*features);
  cout << "Features have been computed for cloud" << endl;
  
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
  \brief  Perform the EKF state update

  Uses a Multiplicative EKF for the attitude - from 	
  Fundamentals of Spacecraft Attitude Determination and Control
  by Markley, F. Landis; Crassidis, John L
  Chapter 6

  \author Benjamin Morrell
  \date 10 May 2018
*/
void NurbSLAM::processStepEKF(float timestep){

  cout << "\nState before process step is:\n" << ekfState << endl;

  switch (processModel){
    case 0: // constant acceleration
      for (int i=0; i < 3; i++){
        // Update positon with constant acceleration model
        ekfState(i) += ekfState(i+3)*timestep + ekfState(i+6)*std::pow(timestep,2.0)*0.5;
        // Update velocity with constant acceleration model
        ekfState(i+3) += ekfState(i+6)*timestep;
      }

      // Angular error in state - constant

      // Update Jacobian with the timestep
      J.setIdentity();
      J(0,3) = timestep;
      J(1,4) = timestep;
      J(2,5) = timestep;
      J(3,6) = timestep;
      J(4,7) = timestep;
      J(5,8) = timestep;
      J(0,6) = 0.5*std::pow(timestep,2.0);
      J(1,7) = 0.5*std::pow(timestep,2.0);
      J(2,8) = 0.5*std::pow(timestep,2.0);
      break;
    case 1: // Constant velocity model
      for (int i=0; i < 3; i++){
        // Update positon with constant velocity model
        ekfState(i) += ekfState(i+3)*timestep;
      }

      // Angular error in state - constant

      // Update Jacobian with the timestep
      J.setIdentity();
      J(0,3) = timestep;
      J(1,4) = timestep;
      J(2,5) = timestep;
      break;
    case 2: // Zero velocity model

      // Angular error in state - constant

      J.setIdentity();

      break;
  }

  

  // Update P
  if (!bKeepPConstant){
    P = J*P*J.transpose() + Q;
  }
  
  cout << "State after process step with timestep " << timestep << "s is:\n" << ekfState << endl;
  cout << "Covariance, P is:\n" << P << endl;

  // Update state - pose
  state.matrix().block(0,3,3,1) = ekfState.block(0,0,3,1);
}

/*! 
  \brief  Perform the SLAM filter updates

  \warning currently a placeholder that just uses the first scan

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::updateSLAMFilter(float timestep){
  // Maybe just loop for each observation
  

  //--------------------------------------------------
  //--------------------------------------------------
  // Get the observation(s)
  //--------------------------------------------------
  //--------------------------------------------------
  // Observation
  transformDelta.matrix() = transformationList[0];// list is of 4x4 matrices
  // Eigen::Affine3f deltaRotOnly(Eigen::AngleAxisf(transformDelta.linear()));
  Eigen::Matrix<float,6,1> observation, predicted;
  predicted.setZero();
  Eigen::Vector3f stateTrans;
  stateTrans(0) = state(0,3);
  stateTrans(1) = state(1,3);
  stateTrans(2) = state(2,3);

  // cout << "Transformation linear is:\n" << transformDelta.linear() << endl;

  // Rotation - get rodrigeuz parameters
  Eigen::Vector3f attitudeErrorRodrigeuz;
  // Axis angle
  Eigen::AngleAxisf attitudeErrorAA(transformDelta.linear());
  // Compute Rodrigeuz parameter
  attitudeErrorRodrigeuz = attitudeErrorAA.axis()*std::tan(attitudeErrorAA.angle()/2.0);

  observation(3) = attitudeErrorRodrigeuz(0);
  observation(4) = attitudeErrorRodrigeuz(1);
  observation(5) = attitudeErrorRodrigeuz(2);
  
  // Translation 
  Eigen::Vector3f deltaTrans; deltaTrans.setZero();
  deltaTrans(0) = transformDelta(0,3);
  deltaTrans(1) = transformDelta(1,3);
  deltaTrans(2) = transformDelta(2,3);

  // Get t_delt = R_delt t_1 + t_delt - t_1
  deltaTrans = deltaTrans + transformDelta.linear()*stateTrans - stateTrans;
  // deltaTrans = deltaTrans + deltaRotOnly*stateTrans - stateTrans;
  
  
  // Put in the observation
  // observation.block(0,0,3,2) = deltaTrans;
  observation(0) = deltaTrans(0);
  observation(1) = deltaTrans(1);
  observation(2) = deltaTrans(2);

  cout << "Observation is:\n" << observation << endl;
  // cout << "Predicted is:\n" << predicted << endl;
  // cout << "Innovation is: \n" << (observation - predicted) << endl;

  //--------------------------------------------------
  //--------------------------------------------------
  // Update R
  //--------------------------------------------------
  //--------------------------------------------------
  
  float angularError = std::abs(attitudeErrorAA.angle());
  float linearError = 0.0;
  for (int k = 0; k < 3; k++){linearError += std::pow(observation(k),2.0);}
  linearError = std::sqrt(linearError);
  cout << "Transform size metrics are, linear: " << linearError << ", angular: " << angularError << endl;

  // Thesholds
  float angularErrorMult = 1.0;
  float linearErrorMult = 1.0;  

  if (angularError > rejectCriteria[0]){
    // High uncertainty if more than 45 degrees
    angularErrorMult = 1e10;
    cout << "Angular error > 45 deg, placing large uncertainty" << endl;
    bRejectAlignment = true;
  }
  if (angularError > rejectCriteria[1]){
    // High uncertainty if more than 90 degrees
    cout << "Angular error > 90 deg, placing very large uncertainty" << endl;
    angularErrorMult = 1e15;
    bRejectAlignment = true;
  }

  if (linearError > rejectCriteria[2]){
    // High uncertainty if more than 5 m NEED TO ADJUST THIS FOR DIFFERENT SCALES
    cout << "Linear error > " << rejectCriteria[2] << " m, placing large uncertainty" << endl;
    linearErrorMult = 1e10;
    bRejectAlignment = true;
  }

  if (linearError > rejectCriteria[3]){
    // High uncertainty if more than 5 m NEED TO ADJUST THIS FOR DIFFERENT SCALES
    cout << "Linear error > " << rejectCriteria[3] << " m, placing very large uncertainty" << endl;
    linearErrorMult = 1e15;
    bRejectAlignment = true;
  }

  if (inlierFractionList[0] < rejectCriteria[4]){
    // Increase uncertainty if fraction is low
    cout << "Very low inlier fraction. Placing large uncertainty" << endl;
    linearErrorMult *= 1e5;
    angularErrorMult *= 1e5;
    bRejectAlignment = true;
  }

  int desiredSize = mp.numRowsDesired*mp.numColsDesired;
  if ((float)numberOfPointsInAlignment/(float)desiredSize < rejectCriteria[5]){
    // High penaty if there are not many points
    cout << "Very low number of points. Placing large uncertainty" << endl;
    linearErrorMult *= 1e5;
    angularErrorMult *= 1e5;
    bRejectAlignment = true;
  }

  if (bRejectAlignment){
    cout << "No update to state, alignment rejected" << endl;
    transformDelta = Eigen::Affine3f::Identity();
    // State not updated
    return;
  }

  linearErrorMult = std::min((float)1e25,linearErrorMult);
  angularErrorMult = std::min((float)1e25,angularErrorMult);

  cout << "Linear Error Mult is: " << linearErrorMult << endl;
  cout << "Angular Error Mult is: " << angularErrorMult << endl;
  
  float noiseObsMultPosSize = noiseObsMultPos*0.1;
  // float noiseObsMultPosErr = noiseObsMultPos*0.5;
  float noiseObsMultAngSize = noiseObsMultAng*0.1;
  // float noiseObsMultAngErr = noiseObsMultAng*0.5;

  // Set R from inlier fraction
  float metric = inlierFractionList[0];
  float threeSig = noiseObsBasePos + noiseObsMultPos*(std::pow(1.0-metric,2.0)) + noiseObsMultPosSize*(std::pow(1.0 - (float)numberOfPointsInAlignment/(float)desiredSize,2.0)) + noiseObsMultPosErr*(std::pow(linearError,4.0));
  R.setIdentity();
  R(0,0) = std::sqrt(threeSig/3.0);
  R(1,1) = std::sqrt(threeSig/3.0);
  R(2,2) = std::sqrt(threeSig/3.0);
  // Angles
  threeSig = noiseObsBaseAng + noiseObsMultAng*(std::pow(1.0-metric,2.0)) + noiseObsMultAngSize*(std::pow(1.0 - (float)numberOfPointsInAlignment/(float)desiredSize,2.0)) + noiseObsMultAngErr*(std::pow(angularError,4.0));
  R(3,3) = std::sqrt(threeSig/3.0);
  R(4,4) = std::sqrt(threeSig/3.0);
  R(5,5) = std::sqrt(threeSig/3.0);

  if (bRejectAlignment){
    R = R*(float)1e25;
  }else{
    R = R*rMatMultiplier*linearErrorMult*angularErrorMult;
  }
  

  cout << "R is: \n" << R << endl;
  // TODO - Revise and tune these settings 


  // Compute the Kalman gain
  K = P*Jh.transpose()*(Jh*P*Jh.transpose() + R).inverse();

  cout << "Kalman Gain is:\n" << K << endl;

  // KF Update 
  Eigen::Matrix<float,12,1> ekfUpdate;

  ekfUpdate = K*(observation - predicted);

  cout << "EKF Update is:\n" << ekfUpdate << endl;

  // Update linear
  ekfState.block(0,0,9,1) += ekfUpdate.block(0,0,9,1);
  // cout << "EKF State is:\n" << ekfState << endl;

  // Output angular
  attitudeErrorRodrigeuz = ekfUpdate.block(9,0,3,1);

  cout << "Rodrigeuz state is:\n" << attitudeErrorRodrigeuz << endl;
  cout << "Norm is: " << 2.0*std::atan(attitudeErrorRodrigeuz.norm()) << "\nVector is: \n" << attitudeErrorRodrigeuz.normalized() << endl;

  cout << "TransformDelta from alignment is:\n" << transformDelta.matrix() << endl;

  // Update state estimate
  cout << "State, pre transformation is:\n" << state.matrix() << endl;  
  
  // Apply rotation 
  Eigen::Affine3f deltaRotTransform(Eigen::AngleAxisf(2.0*std::atan(attitudeErrorRodrigeuz.norm()), attitudeErrorRodrigeuz.normalized()));
  // state.rotate(Eigen::AngleAxisf(2*std:atanf(attitudeErrorRodrigeuz.norm()), attitudeErrorRodrigeuz.normalized()));

  state = deltaRotTransform*state;

  // Apply Translation
  state(0,3) = ekfState(0);
  state(1,3) = ekfState(1);
  state(2,3) = ekfState(2);
  
  cout << "State after transform is:\n" << state.matrix() << endl;

  // Get updated transform delta
  transformDelta = state*oldState.inverse();

  cout << "TransformDelta is:\n" << transformDelta.matrix() << endl;

  // Update covariance
  if (!bKeepPConstant){
    P = (Eigen::Matrix<float, 12,12>::Identity() - K*Jh)*P;
  }
  
  cout << "Covariance P after update step is:\n" << P << endl;

}

/*! 
  \brief  Update the alignment of meshes and update the map

  the "update Map" step

  \author Benjamin Morrell
  \date 30 April 2018
*/
void NurbSLAM::alignAndUpdateMeshes(){
  
  if (bRejectAlignment){
    cout << "Rejected alignment. Not updating the mesh" << endl;
    bRejectAlignment = false;
    return;
  }
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed (new pcl::PointCloud<pcl::PointNormal>(mp.numRowsDesired, mp.numColsDesired));
  
  std::vector<float> searchMetrics;

  Eigen::Affine3f localTransformDelta;
  localTransformDelta = transformDelta;

  int updateID = -1;

  cout << "In Align and update Meshes. TransformDelta EKF is:\n" << transformDelta.matrix() << endl;
  
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
    if (bUseFullAlignmentTransformInUpdate && transformationList.size() > 0){
      cout << "Resetting transform delta for aligning meshes to full alignment" << endl;
      localTransformDelta.matrix() = transformationList[i];
      cout << "In Align and update Meshes. TransformDelta full is:\n" << transformationList[i] << endl;
    }else{
      localTransformDelta = transformDelta;
    }
    
    if (objIDList[i] == -1 && bUseOldStateForNewObjects){
      // New object - use transform before ekf process step
      localTransformDelta = previousState*oldState.inverse(); // Will need to check this...
      cout << "Using old state for add new object. Transform is:\n" << localTransformDelta.matrix() << endl;
    }

    // Align scans with new transformation delta (delat from existing transform)
    pcl::transformPointCloud(*objectMeshList[i], *cloudTransformed, localTransformDelta);

    // Compute normals if not already computed
    if (!bObjectNormalsComputed){
      computeNormals(cloudTransformed);
      bObjectNormalsComputed = true;
    }


    if (objIDList[i] == -1){
      // New object
      cout << "\n\t\t ADDING NEW OBJECT \n\n";
      searchMetrics = mp.computeSearchMetrics(cloudTransformed);
      mp.addObject(cloudTransformed, searchMetrics);

      // ID for the latest object
      updateID = mp.objectMap.size() - 1;

      // std::string filename = "/home/bjm/SpaceCRAFT/ros_ws/updatedSurface_" + std::to_string(updateCount) + ".pcd";
      // mp.writeObjectPCDFile(filename.c_str(), updateID, 125, 125);
      // updateCount++;
    }else{
      cout << "\n\t\t UPDATING OBJECT " << objIDList[i] << "...\n\n";
      // Update object
      // Get size before update
      int nRows = mp.objectMap[objIDList[i]].ctrlPnts().rows();
      int nCols = mp.objectMap[objIDList[i]].ctrlPnts().cols();
      oldObj = mp.objectMap[objIDList[i]];

      mp.updateObject(objIDList[i], cloudTransformed); 

      cout << "Updated object (in nubrSLAM.cpp)" << endl;

      // track if there were any updates - change in the number of control points
      if (nRows != mp.objectMap[objIDList[i]].ctrlPnts().rows() ||
      nCols != mp.objectMap[objIDList[i]].ctrlPnts().cols()){
        cout << "Surface was updated" << endl;
        updateID = objIDList[i];

        // std::string filename = "/home/bjm/SpaceCRAFT/ros_ws/updatedSurface_" + std::to_string(updateCount) + ".pcd";
        // mp.writeObjectPCDFile(filename.c_str(), objIDList[i], 125, 125);
        // updateCount++;

        if (mp.objectMap[objIDList[i]].nansInObject()){
          cout << "Nans in combined object. Reverting" << endl;
          mp.updateObjectInMap(objIDList[i],oldObj);
        }
        mapExtendCount[objIDList[i]] = 0;
      }else{
        // No update
        mapExtendCount[objIDList[i]]++;
        // if (mp.ss.extendDirection[0] != 'N'){
        //   // If an extension computed, but not enough data points
        //   cout << "Incrementing mapExtendCount for object " << objIDList[i] << ", count is now: " << mapExtendCount[objIDList[i]] << endl;
        //   mapExtendCount[objIDList[i]]++;
        // }else{
        //   // No extension computed (e.g. all overlap)
        //   cout << "Map extend count for object " << objIDList[i] << ", at 0" << endl;
        //   mapExtendCount[objIDList[i]] = 0;
        // }
        
        // if (mapExtendCount[objIDList[i]] > mapExtendThreshold){
        //   cout << "Creating new object as Map Extend count for object " << objIDList[i] << " exceeds limit: " << mapExtendThreshold << endl;
        //   // Set ID for new object
        //   updateID = mp.objectMap.size() - 1;
        // }else{
        //   // Flag for no update
        //   updateID = -1;
        // }

        // Flag for no update
          updateID = -1;
      }
      
    }

    if (!bMappingModeOn){
      // Update mesh and feature for the items in the map
      updatePointCloudAndFeaturesInMap(updateID);
    }else{
      cout << "In Mapping mode. data points or map features for localisation to be computed" << endl;
      if (updateID == mp.objectMap.size() - 1){
        // If there is a new object - add to the list
        mapMeshList.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr (new pcl::PointCloud<pcl::PointNormal>));
        mapFeatureList.push_back(pcl::PointCloud<pcl::FPFHSignature33>::Ptr (new pcl::PointCloud<pcl::FPFHSignature33>));
        mapMatchCount.push_back(0);
        mapExtendCount.push_back(0);
      }
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
    mapMatchCount.push_back(0);
    mapExtendCount.push_back(0);
  }else{
    // write a new object over the old
    cout << "update object - replace the mesh and features" << endl;
    mapMeshList[objID] = pcl::PointCloud<pcl::PointNormal>::Ptr (new pcl::PointCloud<pcl::PointNormal>(mtSurf, msSurf, pcl::PointNormal()));
    mapFeatureList[objID] = pcl::PointCloud<pcl::FPFHSignature33>::Ptr (new pcl::PointCloud<pcl::FPFHSignature33>);
  }
  
  // Generate point cloud for NURBS
  mp.pointCloudFromObject3D(objID, msSurf, mtSurf, mapMeshList[objID]);

  if (doesPointCloudHaveNans(mapMeshList[objID])){
    cout << "Nans in point cloud, revert to old" << endl;
    mp.updateObjectInMap(objID,oldObj);
    mp.pointCloudFromObject3D(objID, msSurf, mtSurf, mapMeshList[objID]);
  }

  // Search tree
  // pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_(new pcl::search::KdTree<pcl::PointNormal>);

  cout << "Normal before at (0,0) from NURBS is: " << mapMeshList[objID]->at(0,0).normal_x << " , " << mapMeshList[objID]->at(0,0).normal_y << " , " << mapMeshList[objID]->at(0,0).normal_z << endl;

  // // Estimate normals
  computeNormals(mapMeshList[objID]);
  // // pcl::NormalEstimation<pcl::PointNormal, pcl::PointNormal> nest;
  // // nest.setRadiusSearch(pclNormalRadiusSetting);
  // // // nest.setKSearch(7);
  // // nest.setSearchMethod(search_method_);
  // // nest.setViewPoint(state.matrix()(0,3),state.matrix()(1,3),state.matrix()(2,3)); // Set the viewpoint to the current pose
  // // nest.setInputCloud(mapMeshList[objID]);
  // // nest.compute(*mapMeshList[objID]);
  cout << "computed normals for map cloud" << endl;

  cout << "Normal after at (0,0) from pcl is: " << mapMeshList[objID]->at(0,0).normal_x << " , " << mapMeshList[objID]->at(0,0).normal_y << " , " << mapMeshList[objID]->at(0,0).normal_z << endl;

  // Estimate features
  computeFeatures(mapMeshList[objID], mapMeshList[objID], mapFeatureList[objID]);
  // pcl::FPFHEstimationOMP<pcl::PointNormal,pcl::PointNormal,pcl::FPFHSignature33> fest;
  // fest.setRadiusSearch(pclFeatureRadiusSetting);
  // // fest.setKSearch(7);
  // fest.setInputCloud(mapMeshList[objID]);
  // fest.setSearchMethod(search_method_);
  // fest.setInputNormals(mapMeshList[objID]);
  // fest.compute (*mapFeatureList[objID]);
  cout << "Features have been computed for map cloud" << endl;

  // Set flag that there have been updates
  bMapUpdatedFromScan = true;

}

void NurbSLAM::setInitEKFStates(){
  cout << "Setting up EKF matrices" << endl;
  // EKF init
  P.setIdentity();
  P(0,0) = pNoisePos;
  P(1,1) = pNoisePos;
  P(2,2) = pNoisePos;
  P(3,3) = pNoiseVel;
  P(4,4) = pNoiseVel;
  P(5,5) = pNoiseVel;
  P(6,6) = pNoiseAccel;
  P(7,7) = pNoiseAccel;
  P(8,8) = pNoiseAccel;
  P(9,9) = pNoiseAng;
  P(10,10) = pNoiseAng;
  P(11,11) = pNoiseAng;

  P = P*pNoiseMultiplier;

  Q = P*qNoiseMultiplier;
  // Q.setIdentity();
  // Q(0,0) = pNoisePos;
  // Q(1,1) = pNoisePos;
  // Q(2,2) = pNoisePos;
  // Q(3,3) = pNoiseVel;
  // Q(4,4) = pNoiseVel;
  // Q(5,5) = pNoiseVel;
  // Q(6,6) = pNoiseAccel;
  // Q(7,7) = pNoiseAccel;
  // Q(8,8) = pNoiseAccel;
  // Q(9,9) = pNoiseAng;
  // Q(10,10) = pNoiseAng;
  // Q(11,11) = pNoiseAng;

  J.setIdentity();
  J(0,3) = 1.0;J(1,4) = 1.0;J(2,5) = 1.0;
  J(0,6) = 1.0;J(1,7) = 1.0;J(2,8) = 1.0;
  J(4,6) = 1.0;J(5,7) = 1.0;J(6,8) = 1.0;
  // These will be multiplies by time parameters DeltaT once we have it

  R.setIdentity(); // This will be adjusted once observations are made

  Jh.setZero();
  Jh(0,0) = 1.0;
  Jh(1,1) = 1.0;
  Jh(2,2) = 1.0;
  Jh(3,9) = 1.0;
  Jh(4,10) = 1.0;
  Jh(5,11) = 1.0;

  // init EKF state
  ekfState.setZero();

  cout << "P init is:\n" << P << endl;
}

//-------------------------------------------------------------------
// ------------- Get and set  -----------------------------
//-------------------------------------------------------------------


/*! 
  \brief  Sets the state
  
  Used for initialisation (if not starting at zero)
  Also used in mapping paradigms

  \author Benjamin Morrell
  \date 30 April 2018
*/
  void NurbSLAM::setState(Eigen::Affine3f inputState){
  this->state = inputState;
  oldState = state;
  // Update EKF state
  ekfState(0) = state.matrix()(0,3);
  ekfState(1) = state.matrix()(1,3);
  ekfState(2) = state.matrix()(2,3);

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

bool NurbSLAM::doesPointCloudHaveNans(pcl::PointCloud<pcl::PointNormal>::Ptr cloud){

  for (int i = 0; i < cloud->height; i++){
    for (int j = 0; j < cloud->width; j++){
      if (!pcl::isFinite(cloud->at(j,i))){
        cout << "Nans in point cloud. First at (" << i << ", " << j << ")\n";
        return true;
      }
    }
  }

  return false;
}

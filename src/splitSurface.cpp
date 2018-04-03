#include "cans/splitSurface.h"

using namespace PLib ; 




SplitSurface::SplitSurface(): correspondences_(new pcl::Correspondences), nExtraNew(0), maxDistThreshMultiplier(0.3)
{//corr_filtPtr(new pcl::Correspondences),
  // TDBM May have to set size of dynamic arrays

  newDataIndices.setZero(2,2);
  overlapDataIndices.setZero(2,2);

  nRowColCount.setZero(1,4);
  

}

void SplitSurface::setInputMap(pcl::PointCloud<pcl::PointNormal>::Ptr mapCloud){
  this->mapCloud = mapCloud;
  // TDBM - This may cause issues!
}

void SplitSurface::setInputObservation(pcl::PointCloud<pcl::PointNormal>::Ptr obsCloud){
  this->obsCloud = obsCloud;
  
}

int SplitSurface::getCloudWidth(){
  return mapCloud->width;
}

SplitSurface::~SplitSurface(){
  // Do I need to something here to no leak memory?
  // TDBM
}


/*! 
  \brief  Runs all the steps to determine the split of the new surface and the extension directions

  \param mapCloud  the point cloud for the stored object
  \param obsCloud  the point cloud for the new observation
 
  \author Benjamin Morrell
  \date 2 April 2018
*/
void SplitSurface::splitNewSurfaceObservation(pcl::PointCloud<pcl::PointNormal>::Ptr mapCloud, pcl::PointCloud<pcl::PointNormal>::Ptr obsCloud){

  // Load the surfaces
  setInputMap(mapCloud);
  setInputObservation(obsCloud);

  // Find new (non overlapping) data
  findNonOverlappingData();

  // Get the direction to take the new data
  getNewDataExtendDirection();

  // Get the direction to extend the surface
  getMapDataExtendDirection();

  cout << "Extend Direction is: " << extendDirection << endl;  

}

/*! 
  \brief  Matches clouds, and rejects correspondences to find overlapping and non-overlapping points

  \param mapCloud  the point cloud for the stored object
  \param obsCloud  the point cloud for the new observation
  \param[out] newPointsArray a boolean array with 1's where data is not overlapped in obsCloud
  
  \author Benjamin Morrell
  \date 1 April 2018
*/
// void SplitSurface::findNonOverlappingData(pcl::PointCloud<pcl::PointXYZ>::Ptr mapCloud, pcl::PointCloud<pcl::PointXYZ>::Ptr obsCloud){
void SplitSurface::findNonOverlappingData(){  
  // Estimate normals
  pcl::NormalEstimationOMP<pcl::PointNormal, pcl::PointNormal> nest;
  nest.setRadiusSearch(0.025);
  nest.setInputCloud(obsCloud);
  nest.compute(*obsCloud);
  nest.setInputCloud(mapCloud);
  nest.compute(*mapCloud);

  // Estimate correspondences
  pcl::registration::CorrespondenceEstimationNormalShooting<pcl::PointNormal, pcl::PointNormal, pcl::PointNormal> corrEst;
  corrEst.setInputSource (obsCloud);// Gets a corresponding point for every point in the Source
  corrEst.setSourceNormals (obsCloud);
  corrEst.setInputTarget (mapCloud); // Target for the source to find corresponding points in
  // Test the first 10 correspondences for each point in source, and return the best
  corrEst.setKSearch (10);
  
  corrEst.determineCorrespondences (*correspondences_);

  // Correspondence rejection
  // pcl::Correspondences& corr_filt = *corr_filtPtr;
  pcl::Correspondences& corr_filt = *correspondences_;
  // pcl::registration::CorrespondenceRejectorDistance corrRej;
  // corrRej.setMaximumDistance( 0.05); 
  // corrRej.setInputCorrespondences( correspondences_);
  // corrRej.getCorrespondences( *corr_filtPtr);

  // Create boolean array with 1's where the correspondence was rejected
  newPointsArray.resize(obsCloud->height, obsCloud->width);// Maybe do when load the cloud
  newPointsArray.setZero(obsCloud->height, obsCloud->width);

  float maxDistance = computeDistanceThreshold(); 
  maxDistance *= maxDistThreshMultiplier;
  cout << "Distance threshold is: " << maxDistance << endl;

  int offset;

  for (int i = 0; i < obsCloud->height; i++){
    for (int j = 0; j < obsCloud->width; j++){
      offset = i*obsCloud->width + j; // Row major
      // offset = j*obsCloud->height + i; // Col major
      // cout << "offset is: " << offset << ", query index is: " << corr_filt[offset].index_query << endl;
      if (corr_filt[offset].distance > maxDistance){
        // Not a correspondence, hence is a new point
        // cout << "new point" << endl;
        newPointsArray(i,j) = true;
      }
      // else leave at zero
    }
  }

  cout << "Number of new points: " << newPointsArray.count() << endl;
}

//-------------------------------------------------------------------
/*! 
  \brief  Computes a measure of the mean delta between points, to use as a distance threshold for correspondence

  \param[out] meanDelta - the computed mean delta between adjacent points
  
  \author Benjamin Morrell
  \date 3 April 2018
*/
float SplitSurface::computeDistanceThreshold(){

  Eigen::Vector3f delta;
  float meanDelta = 0;
  int pointCount = 0;

  // Take one row and get the mean distance between points
  for (int j = 0; j < mapCloud->width-2; j++){
    delta = (mapCloud->at(j+1,mapCloud->height/2-1).getArray3fMap() - mapCloud->at(j,mapCloud->height/2-1).getArray3fMap());
    meanDelta += delta.norm();
    pointCount++;
  }
  // Take one column and get the mean distance between points
  for (int i = 0; i < mapCloud->height-2; i++){
    delta = (mapCloud->at(mapCloud->width/2-1,i+1).getArray3fMap()-mapCloud->at(mapCloud->width/2-1,i).getArray3fMap());
    meanDelta += delta.norm();
    pointCount++;
  }

  // Take mean from all the deltas
  meanDelta /= (float)pointCount;

  return meanDelta;
}

//-------------------------------------------------------------------
/*! 
  \brief  Determines the direction for the new data to be take to add to the surface, and the indices for the new data

  \param mapCloud  the point cloud for the stored object
  \param obsCloud  the point cloud for the new observation
  
  \author Benjamin Morrell
  \date 1 April 2018
*/
void SplitSurface::getNewDataExtendDirection(){
  // Get counts of new points in the rows and columns
  newPointsRows.resize(newPointsArray.rows(),1);
  newPointsCols.resize(1,newPointsArray.cols());
  newPointsRows = newPointsArray.rowwise().count().cast<int>();
  newPointsCols = newPointsArray.colwise().count().cast<int>();
  cout << "New Points Rows: " << newPointsRows << endl;
  cout << "New Points Cols: " << newPointsCols << endl;

  // Dimensions of data
  int ms = newPointsArray.rows();
  int mt = newPointsArray.cols();

  // Reset to zero
  nRowColCount.setZero(1,4);
  // COunt will be for L, R, D, U

  getEndNewRowCounts();
  getEndNewColCounts();

  cout << "nRowColCount is [L, R, D, U]: [" << nRowColCount << "]" << endl;

  int extendID;

  int maxCount = nRowColCount.maxCoeff(&extendID);

  cout << "Max ID is: " << extendID << endl;
  cout << "Maximum is: " << maxCount << endl;

  if (maxCount < 2){
    // No new data criteria
    // TODO - TDBM
    cout << "Maximum count is < 2. No extension. All overlap" << endl;

    // No extension, all overlap
    extendDirection = "NN";
    newDataIndices.setZero(2,2);
    overlapDataIndices(0,0) = 0.0;
    overlapDataIndices(1,0) = ms-1;
    overlapDataIndices(0,1) = 0;
    overlapDataIndices(1,1) = mt-1;
    return;
  }

  pcl::Correspondences& corr = *correspondences_;

  if ((nRowColCount==maxCount).count() > 1){
    // Select one
    cout << "Multiple maximums" << endl;
    // TDBM - test this

    // Use the distance from the correspondences, L, R, D, U 
    std::vector<int> idTest(4,0);
    // Test the mid point of each edge
    idTest[0] = mt*ms/2;  // Left
    idTest[1] = mt*ms/2 + mt-1; // Right
    idTest[2] = mt*(ms-1) + mt/2; // Down
    idTest[3] = mt/2; // Up

    float maxD = 0;
    int maxID = -1;

    // Get the distance from the first correspondence search
    for (int i = 1; i < 4; i++){
      if (corr[idTest[i]].distance > maxD){
        // Set to the max if the max is exceeded
        maxD = corr[idTest[i]].distance;
        maxID = i;
      }
    }

    // Take the ID that is the maximum distance away
    extendID = idTest[maxID];
    
  }else{
    cout << "one maximum " << endl;
  }

  

  if (maxCount == ms || maxCount == mt){
    nExtraNew = 0;
  }else{
    nExtraNew = 1;
  }

  cout << "nExtraNew is: " << nExtraNew << endl;

  // Set the indices and search direction. For indices: (0,0) = start in i, (1,0) = end in i, (0,1) is start in j, (1,1) is end in j
  switch (extendID){
    case 0 : // Take the left data
            cout << "Inside Switch statement case 0" << endl;
            // Indices for the new data 
            newDataIndices(0,0) = 0;
            newDataIndices(1,0) = ms - 1;
            newDataIndices(0,1) = 0;
            newDataIndices(1,1) = maxCount + nExtraNew - 1;
            // Indices for the overlapped data
            overlapDataIndices(0,0) = 0;
            overlapDataIndices(1,0) = ms-1;
            overlapDataIndices(0,1) = maxCount + nExtraNew - 1;
            overlapDataIndices(1,1) = mt - 1;
            // Extend direction
            extendDirection[1] = 'L';
            break;
    case 1 : // Take the Right data
            cout << "Inside Switch statement case 1" << endl; 
            // Indices for the new data
            newDataIndices(0,0) = 0;
            newDataIndices(1,0) = ms -1;
            newDataIndices(0,1) = mt - (maxCount + nExtraNew); // take last (maxCount-nExtraNew) columns
            newDataIndices(1,1) = mt - 1;
            // Indices for the overlapped data
            overlapDataIndices(0,0) = 0;
            overlapDataIndices(1,0) = ms-1;
            overlapDataIndices(0,1) = 0;
            overlapDataIndices(1,1) = mt - (maxCount + nExtraNew); 
            // Extend direction
            extendDirection[1] = 'R';
            break;
    case 2 : // Take the Down data
            cout << "Inside Switch statement case 2" << endl;
            // Indices for the new data
            newDataIndices(0,0) = 0;
            newDataIndices(1,0) = maxCount + nExtraNew - 1;
            newDataIndices(0,1) = 0;
            newDataIndices(1,1) = mt - 1;
            // Indices for the overlapped data
            overlapDataIndices(0,0) = maxCount + nExtraNew - 1;
            overlapDataIndices(1,0) = ms - 1;
            overlapDataIndices(0,1) = 0;
            overlapDataIndices(1,1) = mt - 1;
            // Extend direction
            extendDirection[1] = 'D';
            break;
    case 3 : // Take the Up data
            cout << "Inside Switch statement case 3" << endl;
            // Indices for the new data
            newDataIndices(0,0) = ms - (maxCount + nExtraNew);
            newDataIndices(1,0) = ms - 1;
            newDataIndices(0,1) = 0;
            newDataIndices(1,1) = mt - 1;
            // Indices for the overlapped data
            overlapDataIndices(0,0) = 0;
            overlapDataIndices(1,0) = ms - (maxCount + nExtraNew);
            overlapDataIndices(0,1) = 0;
            overlapDataIndices(1,1) = mt - 1;
            // Extend direction
            extendDirection[1] = 'U';
  }

  cout << "Indices out for new are:\n" << newDataIndices << endl;
  cout << "Indices out for overlap are:\n" << overlapDataIndices << endl;
}

//-------------------------------------------------------------------
/*! 
  \brief  Counts the number of rows from start and end that are completely "new" data

  
  \author Benjamin Morrell
  \date 1 April 2018
*/
void SplitSurface::getEndNewRowCounts(){

  bool startFlag = true;
  bool endFlag = true;

  int numberOfRows = newPointsRows.rows();
  int pointsInRow = newPointsCols.cols();

   for (int i = 0; i < numberOfRows; i++){
    if (startFlag){
      if (newPointsRows[i] == pointsInRow){
        // Completely new row
        // Add to count for Down
        nRowColCount[2]++;
      }else{
        startFlag = false;
      }
    }

    if (endFlag){
      if (newPointsRows[numberOfRows - 1 - i] == pointsInRow){
        // completely new row at end
        // Add to count for Up
        nRowColCount[3]++;
      }else{
        endFlag = false;
      }
    }
    if (!startFlag && !endFlag){
      break;
    }
  }
}

//-------------------------------------------------------------------
/*! 
  \brief  Counts the number of Columns from start and end that are completely "new" data

  
  \author Benjamin Morrell
  \date 1 April 2018
*/
void SplitSurface::getEndNewColCounts(){

  bool startFlag = true;
  bool endFlag = true;

  int numberOfCols = newPointsCols.cols();
  int pointsInCol = newPointsRows.rows();
  
  
  for (int i = 0; i < numberOfCols; i++){
    if (startFlag){
      if (newPointsCols[i] == pointsInCol){
        // Completely new Col
        // Add to count for Left
        nRowColCount[0]++;
      }else{
        startFlag = false;
      }
    }

    if (endFlag){
      if (newPointsCols[numberOfCols - 1 - i] == pointsInCol){
        // completely new Col at end
        // Add to count for Right
        nRowColCount[1]++;
      }else{
        endFlag = false;
      }
    }
    if (!startFlag && !endFlag){
      break;
    }
  }
}



//-------------------------------------------------------------------
/*! 
  \brief  Get the direction to extend the surface

  
  \author Benjamin Morrell
  \date 1 April 2018
*/
void SplitSurface::getMapDataExtendDirection(){
  // Dimensions of data
  int ms = newDataIndices(1,0) - newDataIndices(0,0) + 1;
  int mt = newDataIndices(1,1) - newDataIndices(0,1) + 1;

  // int edgeID;

  // std::string directions = "RLUD"; // R switches with L, U switches with D compared to data directions

  // edgeID = directions.find(extendDirection[1]);

  int row;
  int col;

  // Get the opposite edge of the new data
  switch (extendDirection[1]){
    case 'L' : // Get the right edge of the new data
              row = ms/2;
              col = newDataIndices(1,1) - nExtraNew;
              break;
    case 'R' : // Get the Left edge of the new data
              row = ms/2;
              col = newDataIndices(0,1) + nExtraNew;
              break;
    case 'D' : // Get the upper edge of the new data
              row = newDataIndices(1,0) - nExtraNew;          
              col = mt/2;
              break;
    case 'U' : // Get the upper edge of the new data
              row = newDataIndices(0,0) + nExtraNew;          
              col = mt/2;
              break;
    case 'N' : // No extension
              extendDirection = "NN";
              return;
  }
  cout << "Edge indices are (" << row << ", " << col << ")\n";

  // Get the data at that edge
  pcl::PointNormal dataEdgeMid = obsCloud->at(col, row);

  cout << "Point at edge is: " << dataEdgeMid << endl;

  // Surface edges
  int mssrf = mapCloud->height;
  int mtsrf = mapCloud->width;

  // Initialise indice arrays - left is row index, right is column index
  Eigen::Array<int, 4, 2> edgeIndices; // INdices to get edges of surfavce
  Eigen::Array<int, 4, 2> vecPointIndices; // Indices to get just inside the edges to make a vector to the edges

  // L, R, D, U
  // L
  edgeIndices(0,0) = mssrf/2-1;
  edgeIndices(0,1) = 0;
  vecPointIndices(0,0) = mssrf/2-1;
  vecPointIndices(0,1) = 1;
  // R
  edgeIndices(1,0) = mssrf/2-1;
  edgeIndices(1,1) = mtsrf-1;
  vecPointIndices(1,0) = mssrf/2-1;
  vecPointIndices(1,1) = mtsrf-2;
  // D
  edgeIndices(2,0) = 0;
  edgeIndices(2,1) = mtsrf/2-1;
  vecPointIndices(2,0) = 1;
  vecPointIndices(2,1) = mtsrf/2-1;
  // U
  edgeIndices(3,0) = mssrf-1;
  edgeIndices(3,1) = mtsrf/2-1;
  vecPointIndices(3,0) = mssrf-2;
  vecPointIndices(3,1) = mtsrf/2-1;

  Eigen::Vector3f vec1;
  Eigen::Vector3f vec2;

  float theta;
  float dist;
  float minDist = 99999999.9;
  int minID = -1;

  for (int i = 0; i < 4; i ++){

    // Vector from the map cloud edge to the data extend edge
    vec1 = dataEdgeMid.getArray3fMap() - mapCloud->at(edgeIndices(i,1),edgeIndices(i,0)).getArray3fMap();

    // Vector along the map cloud surface to the edge
    vec2 = mapCloud->at(edgeIndices(i,1),edgeIndices(i,0)).getArray3fMap() - mapCloud->at(vecPointIndices(i,1),vecPointIndices(i,0)).getArray3fMap();

    // Angle between the vectors 
    theta = acos(vec2.dot(vec1)/vec1.norm()/vec2.norm());
    cout << "theta for edge " << i << " is: " << theta << endl;

    if (theta > M_PI/2.0f){
      // Do not use the edge
    }else{
      // Compute the distance between the edges
      dist = vec1.norm();

      cout << "Distance to edge " << i << " is: " << dist << endl;

      if (dist < minDist){
        // track the ID for the minimum distance
        minID = i;
        minDist = dist; // Store the ID
      }
    }
  }

  cout << "Minimum distance from edge to edge is: " << minDist << endl;


  // Use the minimum distance ID to select the extension direction for the surface
  switch (minID){
    case 0 : // Left edge is closest, so extend left
            extendDirection[0] = 'L';
            break;
    case 1 : // Right edge is closest, so extend Right
            extendDirection[0] = 'R';
            break;
    case 2 : // Downward edge is closest, so extend down
            extendDirection[0] = 'D';
            break;
    case 3 : // Upper edge is closest, so extend Up
            extendDirection[0] = 'U';
  }
  
  cout << "Map Object Extend Direction is: " << extendDirection[0] << endl;
  

  

  // Loop through each surface edge

}
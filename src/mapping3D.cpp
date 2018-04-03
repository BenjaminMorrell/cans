
#include "cans/mapping3D.h"

using namespace PLib ; 
//-------------------------------------------------------------------
/*! 
  \brief  Default Constructor

  \warning 

  \author Benjamin Morrell
  \date 23 March 2018
*/
Mapping3D::Mapping3D():
    n_ctrl_default(7), order{3,3}, number_downsample_points(45), overlap_d_frac(0.6), 
    search_thresh {0.75,0.75,0.75,0,0,0,0,0,0}, numRowsDesired(45), numColsDesired(45),
    maxNanAllowed(10), removeNanBuffer(0)
{
}



//-------------------------------------------------------------------
// Methods...
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/*! 
  \brief  join two nurbs curve, prioritising the first curves

  \param crv1 the first nurbs curve
  \param crv2 the second nurbs curve
  \param newBeforeOld flag to indicate in curve 2 comes before curve 1
  \param flipKnotParam flag to indicate to flip the knot parametric direction for crv2

  \author Benjamin Morrell
  \date 27 March 2018
*/
NurbsCurvef Mapping3D::joinCurves(NurbsCurvef& crv1, NurbsCurvef& crv2, bool newBeforeOld, bool flipKnotParam){

    // Parameters
    int deg = crv1.degree(); // assumed to be the same for both
    int n_ctrl1 = crv1.ctrlPnts().size();
    int n_ctrl2 = crv2.ctrlPnts().size();
    int n_ctrl_out = n_ctrl1 + n_ctrl2 - deg;
    int n_knot_out = n_ctrl1 + n_ctrl2 + 1;

    // cout << "In Join Curves:\nKnots out: " << n_knot_out << "\tCtrl out: " << n_ctrl_out << endl;

    int j;
    int i;

    // Initialise output knot vector and control points
    Vector_FLOAT knots_out(n_knot_out);
    Vector_HPoint3Df ctrl_out(n_ctrl_out);

    if (newBeforeOld){
        // crv 2 comes first

        if (flipKnotParam){
            // Take crv 2 from the end and reverse direction
            // Knots
            for (i = 0; i < n_knot_out; i++){
                if (i < n_ctrl2){
                    // Take curve 2 knots from the end until the first zero
                    knots_out[i] = (1 - crv2.knot()[n_ctrl2 + deg - i])/2.0;
                }else{
                    // Take curve 1 knots
                    j = deg + i - n_ctrl2;
                    knots_out[i] = crv1.knot()[j]/2.0 + 0.5;
                }
            }
            // Control points
            for (i = 0; i < n_ctrl_out; i++){
                if (i < n_ctrl2 - deg ){
                    // Take curve 2 control points from the end, missing the last #deg control points
                    ctrl_out[i] = crv2.ctrlPnts()[n_ctrl2 - 1 - i];
                }else{
                    // Take curve 1 control points
                    j = i - (n_ctrl2 - deg);
                    ctrl_out[i] = crv1.ctrlPnts()[j]; 
                }
            }
        }else{
            // Take crv2 from the start
            // Knots
            for (i = 0; i < n_knot_out; i++){
                if (i < n_ctrl2 ){
                    // Take curve 2 knots from the start up to the last non-one term
                    knots_out[i] = crv2.knot()[i]/2.0;
                }else{
                    // Take curve 1 knots from the last zero term
                    j = deg + i - (n_ctrl2);
                    knots_out[i] = crv1.knot()[j]/2.0 + 0.5;
                }
            }
            // Control points
            for (i = 0; i < n_ctrl_out; i++){
                if (i < n_ctrl2 - deg ){
                    // Take curve 2 control points missing the last #deg 
                    ctrl_out[i] = crv2.ctrlPnts()[i];
                }else{
                    // Take curve 1 control points
                    j = i - (n_ctrl2 - deg);
                    ctrl_out[i] = crv1.ctrlPnts()[j]; 
                }
            }
        }

    }else{
        // crv 1 comes first

        if (flipKnotParam){
            // Take crv 2 from the start and reverse direction
            // Knots
            for (i = 0; i < n_knot_out; i++){
                if (i < n_ctrl1 + 1){
                    // Take curve 1 knots from the start
                    knots_out[i] = crv1.knot()[i]/2.0;
                }else{
                    // Take curve 2 knots in reverse from the first non-one term
                    j = i - (n_ctrl1+1);
                    knots_out[i] = (1 - crv2.knot()[n_ctrl2 - 1 - j])/2.0 + 0.5;
                }
            }
            // Control points
            for (i = 0; i < n_ctrl_out; i++){
                if (i < n_ctrl1 ){
                    // Take curve 1 control points from the start
                    ctrl_out[i] = crv1.ctrlPnts()[i];
                }else{
                    // Take curve 2 control points missing the last #deg control points
                    j = deg + i - n_ctrl1;
                    ctrl_out[i] = crv2.ctrlPnts()[n_ctrl2 - 1 - j]; 
                }
            }
        }else{
            // Take crv2 from the end
            // Knots
            for (i = 0; i < n_knot_out; i++){
                if (i < n_ctrl1 + 1){
                    // Take curve 1 knots from the start
                    knots_out[i] = crv1.knot()[i]/2.0;
                }else{
                    // Take curve 2 knots from the first non-zero term
                    j = deg + 1 + i - (n_ctrl1+1);
                    knots_out[i] = crv2.knot()[j]/2.0 + 0.5;
                }
            }
            // Control points
            for (i = 0; i < n_ctrl_out; i++){
                if (i < n_ctrl1 ){
                    // Take curve 1 control points from the start
                    ctrl_out[i] = crv1.ctrlPnts()[i];
                }else{
                    // Take curve 2 control points missing the first #deg control points
                    j = deg + i - n_ctrl1;
                    ctrl_out[i] = crv2.ctrlPnts()[j]; // TODO is it inefficient to keep calling this function?
                }
            }

        }
    }

    // Create NURBS Curve
    NurbsCurvef crv_out(ctrl_out,knots_out,deg);

    return crv_out;

}


//-------------------------------------------------------------------
/*! 
  \brief  join two nurbs surfaces prioritising the first surface

  \param srf1 the first nurbs surface
  \param srf2 the second nurbs surface
  \param extendDirection string to indicate how the surfaces are to be joined
  Options are L, R, U, D. Is a string of two values, one for srf1 and one for srf2

  \warning  when parameters are flipped from the second surfce, this is assumed to be from a rotation
    and not from a mirroring
  \warning  Assumes input surfaces have equal numbers of contorl points across dimension perpindicular to how they are joining
  \warning  Assumed the data is organised so that going down a column in the matrix is going "up" the surface parametrically 
  
  \author Benjamin Morrell
  \date 28 March 2018
*/
NurbsSurface<float,3> Mapping3D::joinSurfaces(NurbsSurfacef& srf1, NurbsSurfacef& srf2, char * extendDir){

    int deg = srf1.degreeU();// assume the same degree in each direction
    int row_not_col;
    int n_ctrl1;
    int n_ctrl2;
    bool newBeforeOld;
    bool flipKnotParam;
    Vector_FLOAT knots1;
    NurbsCurvef crv1;
    NurbsCurvef crv2;
    NurbsCurvef crvOut;
    Vector_HPoint3Df ctrlVec1; // TODO maybe don't make copies here - just use the data directly where possible
    Vector_HPoint3Df ctrlVec2;

    Matrix_HPoint3Df ctrlOut;
    Vector_FLOAT knotsOut;

    NurbsSurfacef* srfOut;

    // Get whether srf 2 will be using rows (1) or columns (0)
    if (extendDir[1] == 'L' || extendDir[1] == 'R'){
        row_not_col = true; // use rows
        n_ctrl2 = srf2.ctrlPnts().cols();
    }else{
        row_not_col = false; // use columns
        n_ctrl2 = srf2.ctrlPnts().rows();
    }

    // Different cases for rows or columns for srf1
    if (extendDir[0] == 'L' || extendDir[0] == 'R'){
        // -----------------------------------------------------------------
        // Row-wise expansion
        // -----------------------------------------------------------------
        n_ctrl1 = srf1.ctrlPnts().cols();

        // Init
        ctrlOut.resize(n_ctrl1,n_ctrl1 + n_ctrl2 - deg);
        knotsOut.resize(n_ctrl1 + n_ctrl2 + 1);

        
        knots1 = srf1.knotV();        

        // Determine order of combination
        if (extendDir[0] == 'L'){
            newBeforeOld = true;

            // Determine whether to flip parameters
            if (extendDir[1] == 'L' || extendDir[1] == 'D'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }        
        }else{// Extend Right
            newBeforeOld = false;

            // Determine whether to flip parameters
            if (extendDir[1] == 'R' || extendDir[1] == 'U'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }

        }

        // Loop for each row of control points 
        for (int i = 0; i < n_ctrl1; i++){
            ctrlVec1 = getMatRow(srf1.ctrlPnts(),i);
            crv1.reset(ctrlVec1,knots1,deg);

            if (row_not_col){
                // Use rows
                
                if (flipKnotParam){
                    // take rows in different order - from the end first 
                    ctrlVec2 = getMatRow(srf2.ctrlPnts(),srf2.ctrlPnts().rows() - 1 - i);    
                }else{
                    ctrlVec2 = getMatRow(srf2.ctrlPnts(),i);
                }

                // Set curve 2
                crv2.reset(ctrlVec2,srf2.knotV(),deg);

            }else{
                // Use columns
                
                if (flipKnotParam){
                    // take columns in same order: assume a rotation, so order changes if not flipped
                    ctrlVec2 = getMatCol(srf2.ctrlPnts(),i,false); 
                }else{
                    // take columns in reverse order: assume a rotation, so order changes if not flipped
                    ctrlVec2 = getMatCol(srf2.ctrlPnts(),srf2.ctrlPnts().cols() - 1 - i,false);
                }

                // Set Curve 2
                crv2.reset(ctrlVec2,srf2.knotU(),deg);
            }

            // Join the curves
            crvOut = joinCurves(crv1,crv2,newBeforeOld,flipKnotParam);

            insertMatRow(ctrlOut,crvOut.ctrlPnts(),i);

            knotsOut += crvOut.knot();
        }

        // Average knot vector
        for (int i = 0; i < knotsOut.size(); i++){
            knotsOut[i] /= n_ctrl1;
        }
        
        // create surface
        srfOut = new NurbsSurfacef(deg,deg,knots1,knotsOut,ctrlOut);
        
    }else{
        // -----------------------------------------------------------------
        // Column-wise expansion
        // -----------------------------------------------------------------
        n_ctrl1 = srf1.ctrlPnts().rows();

        // Init
        ctrlOut.resize(n_ctrl1 + n_ctrl2 - deg, n_ctrl1);
        knotsOut.resize(n_ctrl1 + n_ctrl2 + 1);

        knots1 = srf1.knotU();        

        // Determine order of combination
        if (extendDir[0] == 'D'){
            newBeforeOld = true;

            // Determine whether to flip parameters
            if (extendDir[1] == 'D' || extendDir[1] == 'L'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }        
        }else{// Extend Up
            newBeforeOld = false;

            // Determine whether // take columns in same order: assume a rotati// take columns in same order: assume a rotation, so order changes if not flippedon, so order changes if not flippedto flip parameters
            if (extendDir[1] == 'U' || extendDir[1] == 'R'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }

        }

        // Loop for each column of control points 
        for (int j = 0; j < n_ctrl1; j++){
            ctrlVec1 = getMatCol(srf1.ctrlPnts(),j,false);
            crv1.reset(ctrlVec1,knots1,deg);

            if (row_not_col){
                // Use rows
                // TODO: Check if fliplr is needed here...
                if (flipKnotParam){
                    // Assume rotated - so when flipped, the orders align 
                    ctrlVec2 = getMatRow(srf2.ctrlPnts(),j);    
                }else{
                    // Assume rotated - so when not flipped, the order needs to change
                    ctrlVec2 = getMatRow(srf2.ctrlPnts(),srf2.ctrlPnts().rows() - 1 - j);
                }
                
                // Set cruve 2
                crv2.reset(ctrlVec2,srf2.knotV(),deg);

            }else{
                // Use columns
                
                if (flipKnotParam){
                    // take columns in different order - from the end first 
                    ctrlVec2 = getMatCol(srf2.ctrlPnts(),srf2.ctrlPnts().cols() - 1 - j, false); 
                }else{
                    ctrlVec2 = getMatCol(srf2.ctrlPnts(),j, false);
                }

                // Set curve 2 
                crv2.reset(ctrlVec2,srf2.knotU(),deg);
            }

            
            // Join the curves,true,true
            crvOut = joinCurves(crv1,crv2,newBeforeOld,flipKnotParam);

            insertMatCol(ctrlOut,crvOut.ctrlPnts(),j, false);

            knotsOut += crvOut.knot();
        }

        // Average knot vector
        for (int i = 0; i < knotsOut.size(); i++){
            knotsOut[i] /= n_ctrl1;
        }
        
        // create surface
        srfOut = new NurbsSurfacef(deg,deg,knotsOut,knots1,ctrlOut);

    }
    

    return *srfOut;

}

//-------------------------------------------------------------------
// ------------- Scan Processing  -----------------------------
//-------------------------------------------------------------------

//-------------------------------------------------------------------

/*! 
  \brief  Searches through an organised pointcloud and makes a boolean array with 1s where the point cloud is nan

  \param nanArray   the output boolean array
  \param cloud      a pointer to the point cloud to process
  
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::meshFromScan(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOut, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudIn){
  // Initialise
  Eigen::Array<bool, Eigen::Dynamic, 1> rowFlags(cloudIn->height, 1); rowFlags.setOnes(cloudIn->height, 1);
  Eigen::Array<bool, 1, Eigen::Dynamic> colFlags(1, cloudIn->width); colFlags.setOnes(1, cloudIn->width);
  // Initialise Nan array
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> nanArray(cloudIn->height,cloudIn->width);
  nanArray.setZero(cloudIn->height,cloudIn->width);

  std::cerr << "PointCloud before filtering: W: " << cloudIn->width << "\tH: " << cloudIn->height << "\tdata points." << std::endl;

  // Fill the Nan Array
  getNanMatrixFromPointCloud(nanArray, cloudIn);

  cout << "Number of Nans: " << nanArray.count() << endl;

  bool exitFlag = false;

  bool nansPresent = true;

  // Loop to remove Nans
  while (!exitFlag){
    nansPresent = removeRowsNan(nanArray, rowFlags);

    if (nansPresent){
      nansPresent = removeColsNan(nanArray, colFlags);
      if (nansPresent){
        // Check dimensions 
        cout << "Nans removed: Row count: " << rowFlags.count() << "\nCol Count: " << colFlags.count() << endl;
        if (rowFlags.count() <= this->numRowsDesired || colFlags.count() <= this->numColsDesired){
          exitFlag = true;
          cout << "Exiting because dimensions are too small" << endl;
        }
      }else{
        exitFlag = true;
        cout << "Exiting because there are few enough NaNs" << endl;
      }
    }else{
      exitFlag = true;
      cout << "Exiting because there are few enough NaNs" << endl;
    }
  }

  // Downsample 
  downsampleRow(rowFlags);
  downsampleCol(colFlags);

  cout << "Number of nans left: " << nanArray.count() << endl;

  // Extract point cloud
  int ii = 0;
  int jj = 0;
  int ijk = 0;
  Eigen::Array<int, 2, Eigen::Dynamic> nanIndices(2,nanArray.count());
  nanIndices.setConstant(2,nanArray.count(),-1);

  for (int i = 0; i < cloudIn->height; i++){
    // Reset jj index
    jj = 0; 
    // For each row selected in rowFlags
    if (rowFlags(i,0)){
      for (int j = 0; j < cloudIn->width; j++){
        // For each column selected in colFlags
        if (colFlags(0,j)){
          // Get the data from the cloud
          cloudOut->at(jj,ii) = cloudIn->at(j,i);

          // Store indices if the value is nan
          if (!pcl::isFinite(cloudIn->at(j,i))){
            nanIndices(0,ijk) = ii;
            nanIndices(1,ijk) = jj;
            ijk++; 
          }

          // Increment the new column index
          jj++;
        }
      }
      // Increment the new row index
      ii++;
    }
  }

  cout << "Finished copying data across" << endl;

  cout << "NanIndices are: " << nanIndices << endl;

  // Average Nans
  averageOutNans(cloudOut, nanIndices);
  
}

/*! 
  \brief  Searches through an organised pointcloud and makes a boolean array with 1s where the point cloud is nan

  \param nanArray   the output boolean array
  \param cloud      a pointer to the point cloud to process
  
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::getNanMatrixFromPointCloud(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){

   // Loop through Point cloud 
  for (int i = 0; i < cloud->height; i++){
    for (int j = 0; j < cloud->width; j++){
      // If nan value
      if (!pcl::isFinite(cloud->at(j,i))){ // cloud->at(col,row)
        nanArray(i,j) = true;
      }else{
        nanArray(i,j) = false;
      }
    }
  }
}

/*! 
  \brief  Removes rows as being actively selected, and then clears out nanArray to no longer consider that row

  \param nanArray       input boolean array with true where there were nan values in a point cloud
  \param rowFlags       the array of boolean flags indicating the active rows
  \param maxNanAllowed  maximum number of NaNs that is permissible

  The function will return true if rows were removed and false if they were not.
  Uses removeNanBuffer to select more rows in one go. 
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
bool Mapping3D::removeRowsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags){

  Eigen::Array<int, Eigen::Dynamic, 1> rowNan = nanArray.rowwise().count().cast<int>();

  // Maximum number of Nans in a row
  int maxNan = rowNan.maxCoeff();

  cout << "Max NaN count: " << maxNan << "\tmax allowed is: " << this->maxNanAllowed << endl;

  // If the maximum number Nans is above the set limit
  if (maxNan > this->maxNanAllowed){
    for (int i = 0; i < rowNan.rows(); i++){
      // If this rows has equal to the maximum number of nans
      if (rowNan(i) >= maxNan-this->removeNanBuffer){
        // Set the row values to zero in the array
        nanArray.row(i) = Eigen::Array<bool, 1, Eigen::Dynamic>::Zero(1,nanArray.cols());

        // Set the rowFlags value to zero
        rowFlags(i) = false;
      }
    }
    return true;
  }else{
    return false;
  }
}

/*! 
  \brief  Removes cols as being actively selected, and then clears out nanArray to no longer consider that col

  \param nanArray       input boolean array with true where there were nan values in a point cloud
  \param colFlags       the array of boolean flags indicating the active cols
  \param maxNanAllowed  maximum number of NaNs that is permissible

  The function will return true if cols were removed and false if they were not.
  Uses removeNanBuffer to select more cols in one go. 
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
bool Mapping3D::removeColsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags){

  Eigen::Array<int, 1, Eigen::Dynamic> colNan = nanArray.colwise().count().cast<int>();

  // Maximum number of Nans in a col
  int maxNan = colNan.maxCoeff();

  cout << "Max NaN count: " << maxNan << "\tmax allowed is: " << this->maxNanAllowed << endl;

  // If the maximum number Nans is above the set limit
  if (maxNan > this->maxNanAllowed){
    for (int i = 0; i < colNan.cols(); i++){
      // If this cols has equal to the maximum number of nans
      if (colNan(i) >= maxNan-this->removeNanBuffer){
        // Set the col values to zero in the array
        nanArray.col(i) = Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(nanArray.rows(), 1);

        // Set the colFlags value to zero
        colFlags(i) = false;

      }
    }
    return true;
  }else{
    return false;
  }
}

/*! 
  \brief  select a set number of rows - indicated by an array of boolean flags

  \param rowFlags       the array of boolean flags indicating the active rows
  \param numRowsDesired The number of rows desired out 
  
  Uses linspace to exclude rows so that the number of 1's in rowFlags is equal to numRowsDesired
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::downsampleRow(Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags){

  // Create linspaced array 
  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArrayf.setLinSpaced(this->numRowsDesired, 0, rowFlags.count()-1);
  selectArray = selectArrayf.round().cast<int>();

  int j = 0;

  // Make the flags false if not included
  for (int i = 0; i < rowFlags.rows(); i++){
    if (rowFlags(i)){
      if (!(selectArray == j).any()){
        rowFlags(i) = false;
      }
      j++;
    }
  }
}

/*! 
  \brief  select a set number of columns - indicated by an array of boolean flags

  \param colFlags       the array of boolean flags indicating the active columns
  \param numColDesired The number of columns desired out 
  
  Uses linspace to exclude columnss so that the number of 1's in colFlags is equal to numColsDesired
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::downsampleCol(Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags){

  // Create linspaced array 
  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArrayf.setLinSpaced(this->numColsDesired, 0, colFlags.count()-1);
  selectArray = selectArrayf.round().cast<int>();
  
  int j = 0;

  // Make the flags false if not included
  for (int i = 0; i < colFlags.cols(); i++){
    if (colFlags(i)){
      if (!(selectArray == j).any()){
        colFlags(i) = false;
      }
      j++;
    }
  }
}

/*! 
  \brief  Steps through nans in a point cloud and replaces with averages from neighbouring points

  \param cloud       The point cloud to modify
  \param nanIndices  An array of indices of the points to replace
  
  Uses linspace to exclude columnss so that the number of 1's in colFlags is equal to numColsDesired
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::averageOutNans(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic> nanIndices){
  cout << "Number of cols: " << nanIndices.cols() << endl;
  for (int i = 0; i < nanIndices.cols(); i++){
    if (nanIndices(0,i) == -1){
      // Flag meaning that those terms are not valid (not used)
      break;
    }
    cout << "Before: " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
    regionAverage(cloud,nanIndices(0,i),nanIndices(1,i));
    cout << "After:  " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
  }

}

/*! 
  \brief  Replaces a given point by projecting from nearby points

  \param cloud  pointer to cloud to process
  \param i      row index of point to change
  \param j      col index of point to change
  
  Will search +_ in i and j, then try to take the delta from i-2 to i-1 and apply to i-1 to get a new point for i 
  (similarly for i+1, j-1, j+1). The results are averaged together. 
  If the neighbouring points are nan, they are not used. 
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::regionAverage(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, int i, int j){

  pcl::PointXYZ average(0.0,0.0,0.0);

  int add_count = 0; // to track how many values are averaged. 

  // Back in i
  if (i > 0){
    if (pcl::isFinite(cloud->at(j,i-1))){
      if (i > 1){
        if (pcl::isFinite(cloud->at(j,i-2))){
          // Step with same delta as neighbours
          average.getArray3fMap() += 2*cloud->at(j,i-1).getArray3fMap() - cloud->at(j,i-2).getArray3fMap();
        }else{
          average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
        }

      }else{
        average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Forward in i
  if (i < cloud->height-1){
    if (pcl::isFinite(cloud->at(j,i+1))){
      if (i < cloud->height-2){
        if (pcl::isFinite(cloud->at(j,i+2))){
          // Step with same delta as neighbours
          average.getArray3fMap() += 2*cloud->at(j,i+1).getArray3fMap() - cloud->at(j,i+2).getArray3fMap();    
        }else{
          average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
        }
      }else{
        average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Back in j
  if (j > 0){
    if (pcl::isFinite(cloud->at(j-1,i))){
      if (j > 1){
        if (pcl::isFinite(cloud->at(j-2,i))){
          // Step with same delta as neighbours
          average.getArray3fMap() += 2*cloud->at(j-1,i).getArray3fMap() - cloud->at(j-2,i).getArray3fMap();
        }else{
          average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
        }
      }else{
        average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Forward in j
  if (j < cloud->width-1){
    if (pcl::isFinite(cloud->at(j+1,i))){
      if (j < cloud->width-2){
        if (pcl::isFinite(cloud->at(j+2,i))){
          // Step with same delta as neighbours
          average.getArray3fMap() += 2*cloud->at(j+1,i).getArray3fMap() - cloud->at(j+2,i).getArray3fMap();
        }else{
          average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
        }
      }else{
        average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Divide by count to get the average
  average.getArray3fMap() /= (float)add_count;

  //update cloud
  cloud->at(j,i) = average;

}




//-------------------------------------------------------------------
// ------------- CONVENIENCE FUNCTIONS -----------------------------
//-------------------------------------------------------------------


//-------------------------------------------------------------------
/*! 
  \brief  Extracts a row of a matrix

  \param data   from which to extract
  \param row_id row to extract from 
  
  \author Benjamin Morrell
  \date 28 March 2018
*/
Vector_HPoint3Df Mapping3D::getMatRow(Matrix_HPoint3Df data,int row_id){
  
  int n_ctrl = data.rows();
  
  Vector_HPoint3Df ctrlRow(n_ctrl);
  
  for (int j = 0; j < n_ctrl; j++){
    ctrlRow[j] = data(row_id,j);
  }

  return ctrlRow;
}


//-------------------------------------------------------------------
/*! 
  \brief  Extracts a column of a matrix

  \param data   from which to extract
  \param col_id column to extract from 
  \param rev    flag to indicate to reverse the direction of sampling
 
  \author Benjamin Morrell
  \date 28 March 2018
*/
Vector_HPoint3Df Mapping3D::getMatCol(Matrix_HPoint3Df data, int col_id, bool rev){
  
  int n_ctrl = data.rows();
  
  Vector_HPoint3Df ctrlRow(n_ctrl);
  
  for (int i = 0; i < n_ctrl; i++){
    if (rev){
        // Take in reverse, as parametric direction goes from bottom to top
        ctrlRow[i] = data(n_ctrl - 1 - i,col_id);
    }else{
        ctrlRow[i] = data(i,col_id);
    }
      
  }

  return ctrlRow;
}

//-------------------------------------------------------------------
/*! 
  \brief  Insert data into row of a matrix

  \param data   matrix to insert row
  \param row    vector of data to insert
  \param row_id row to add data to
  
  \author Benjamin Morrell
  \date 28 March 2018
*/
void Mapping3D::insertMatRow(Matrix_HPoint3Df& data,Vector_HPoint3Df row, int row_id){
  // Function to insert data into an exisitng matrix - across a whole row
  for (int j = 0; j < data.cols(); j++){
    data(row_id,j) = row[j];
  }

}

//-------------------------------------------------------------------
/*! 
  \brief  Insert data into column of a matrix

  \param data   matrix to insert column
  \param col    vector of data to insert
  \param col_id column to add data to
  \param rev    flag to insert in reverse order - for NURBS surfaces
  
  \author Benjamin Morrell
  \date 28 March 2018
*/
void Mapping3D::insertMatCol(Matrix_HPoint3Df& data,Vector_HPoint3Df col, int col_id, bool rev){
  // Function to insert data into an exisitng matrix - across a whole row
  for (int i = 0; i < data.rows(); i++){
      if (rev) {
          data(i,col_id) = col[data.rows() - 1 - i];
      }else{
          data(i,col_id) = col[i];
      }
  }

}
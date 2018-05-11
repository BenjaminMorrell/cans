
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
    nCtrlDefault{15,15}, order{3,3},
    searchThresh(7), numRowsDesired(45), numColsDesired(45),
    maxNanAllowed(10), removeNanBuffer(0), numberOfMetrics(7), msSurf(125), mtSurf(125),
    knotInsertionFlag(true), numInsert(3), deltaKnotInsert(1e-2), newRowColBuffer(0), useNonRectData(false),
    bFilterZ(false), nPointsZLim(400), zThreshMultiplier(0.03), bRejectScan(false),
    minRowsColsAllowed(15), maxNanPercentage(0.1)
{
  searchThresh[0] = 7.75;
  searchThresh[1] = 7.75;
  searchThresh[2] = 7.75;
  for (int i = 3; i < 7; i++){
    searchThresh[i] = 0.0;  
  }
}

//-------------------------------------------------------------------
// Methods...
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// ------------- HGHER LEVEL FUNCTIONALITY  -----------------------------
//-------------------------------------------------------------------

int Mapping3D::processScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Affine3f transform){

  // Initialise data
  // clouds
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudReduced (new pcl::PointCloud<pcl::PointNormal>(numRowsDesired, numColsDesired));
  pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed (new pcl::PointCloud<pcl::PointNormal>(numRowsDesired, numColsDesired));
  // Search parameters
  std::vector<float> searchMetrics;
  int objID; 
  bRejectScan = false;

  // cout << "cloud in at (0,0): " << cloud->at(0,0) << endl;

  // Process Scan 1
  meshFromScan(cloudReduced, cloud);

  if (bRejectScan){
    cout << "Exiting without processsing" << endl;
    return -1;
  }

  // cout << "cloudReduced at (0,0): " << cloudReduced->at(0,0) << endl;

  cout << "transform is: " << transform.matrix() << endl;

  // Transform
  pcl::transformPointCloud(*cloudReduced, *cloudTransformed, transform);

  // cout << "cloud transformed at (0,0): " << cloudTransformed->at(0,0) << endl;

  // Compute Metrics
  searchMetrics = computeSearchMetrics(cloudTransformed);

  // cout << "After compute metrics, cloud transformed at (0,0): " << cloudTransformed->at(0,0) << endl;

  // Data Association
  objID = dataAssociation(searchMetrics);

  // Add new object if there is no match
  if (objID == -1){
    // New object
    addObject(cloudTransformed, searchMetrics);
    return 2; // 2 indicates a new object has been added
  }

  // Update object
  updateObject(objID, cloudTransformed); 


  return 1; // Object updated
}


//-------------------------------------------------------------------
// ------------- SURFACE JOINING  -----------------------------
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

    float factor1 = (float)n_ctrl1/(float)(n_ctrl1 + n_ctrl2);
    float factor2 = (float)n_ctrl2/(float)(n_ctrl1 + n_ctrl2);

    // cout << "In join curves, factor 1 is: " << factor1 << ", factor2 is: " << factor2 << endl;


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
                    knots_out[i] = (1 - crv2.knot()[n_ctrl2 + deg - i])*factor2;
                }else{
                    // Take curve 1 knots
                    j = deg + i - n_ctrl2;
                    knots_out[i] = crv1.knot()[j]*factor1 + factor2;
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
                    knots_out[i] = crv2.knot()[i]*factor2;
                }else{
                    // Take curve 1 knots from the last zero term
                    j = deg + i - (n_ctrl2);
                    knots_out[i] = crv1.knot()[j]*factor1 + factor2;
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
                    knots_out[i] = crv1.knot()[i]*factor1;
                }else{
                    // Take curve 2 knots in reverse from the first non-one term
                    j = i - (n_ctrl1+1);
                    knots_out[i] = (1 - crv2.knot()[n_ctrl2 - 1 - j])*factor2 + factor1;
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
                    knots_out[i] = crv1.knot()[i]*factor1;
                }else{
                    // Take curve 2 knots from the first non-zero term
                    j = deg + 1 + i - (n_ctrl1+1);
                    knots_out[i] = crv2.knot()[j]*factor2 + factor1;
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

    // cout << "Output knots are: " << knots_out << endl;

    // Create NURBS Curve
    NurbsCurvef crv_out(ctrl_out,knots_out,deg);

    return crv_out;

}

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
NurbsSurface<float,3> Mapping3D::joinSurfaces(NurbsSurfacef& srf1, NurbsSurfacef& srf2, std::string extendDirection){
    // TDBM make this templated
    
    int deg = srf1.degreeU();// assume the same degree in each direction
    int row_not_col;
    int n_ctrl1;
    int n_ctrl2;
    int n_ctrlCheck2;
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

    if (srf1.degreeU() != srf2.degreeU() || srf1.degreeV() != srf2.degreeV()){
          cout << "Error: need matching degrees for join direction. Exiting with empty surface" << endl;
          srfOut = new NurbsSurfacef();
          return *srfOut;// This may break
    }

    // Get whether srf 2 will be using rows (1) or columns (0)
    if (extendDirection[1] == 'L' || extendDirection[1] == 'R'){
        row_not_col = true; // use rows
        n_ctrl2 = srf2.ctrlPnts().cols();       
        n_ctrlCheck2 = srf2.ctrlPnts().rows();  
    }else{
        row_not_col = false; // use columns
        n_ctrl2 = srf2.ctrlPnts().rows();
        n_ctrlCheck2 = srf2.ctrlPnts().cols();  
    }

    // Different cases for rows or columns for srf1
    if (extendDirection[0] == 'L' || extendDirection[0] == 'R'){
        // -----------------------------------------------------------------
        // Row-wise expansion
        // -----------------------------------------------------------------
        // Input checks
        if (n_ctrlCheck2 != srf1.ctrlPnts().rows()){
          cout << "Error!: Need to have matched numbers of control points along the join direction\nReturning an empty surface" << endl;
          cout << "Srf1 dimensions are: (" << srf1.ctrlPnts().rows() << ", " << srf1.ctrlPnts().cols() << ")\n";
          cout << "Srf2 dimensions are: (" << srf2.ctrlPnts().rows() << ", " << srf2.ctrlPnts().cols() << ")\n";
          srfOut = new NurbsSurfacef();
          return *srfOut;
        }

        n_ctrl1 = srf1.ctrlPnts().cols();

        // Init
        ctrlOut.resize(n_ctrl1,n_ctrl1 + n_ctrl2 - deg);
        knotsOut.resize(n_ctrl1 + n_ctrl2 + 1);

        
        knots1 = srf1.knotV();        

        // Determine order of combination
        if (extendDirection[0] == 'L'){
            newBeforeOld = true;

            // Determine whether to flip parameters
            if (extendDirection[1] == 'L' || extendDirection[1] == 'D'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }        
        }else{// Extend Right
            newBeforeOld = false;

            // Determine whether to flip parameters
            if (extendDirection[1] == 'R' || extendDirection[1] == 'U'){
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
        // Input checks
        if (n_ctrlCheck2 != srf1.ctrlPnts().cols()){
          cout << "Error!: Need to have matched numbers of control points along the join direction\nReturning an empty surface" << endl;
          srfOut = new NurbsSurfacef();
          return *srfOut;
        }

        n_ctrl1 = srf1.ctrlPnts().rows();

        // Init
        ctrlOut.resize(n_ctrl1 + n_ctrl2 - deg, n_ctrl1);
        knotsOut.resize(n_ctrl1 + n_ctrl2 + 1);

        knots1 = srf1.knotU();        

        // Determine order of combination
        if (extendDirection[0] == 'D'){
            newBeforeOld = true;

            // Determine whether to flip parameters
            if (extendDirection[1] == 'D' || extendDirection[1] == 'L'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }        
        }else{// Extend Up
            newBeforeOld = false;

            // Determine whether // take columns in same order: assume a rotati// take columns in same order: assume a rotation, so order changes if not flippedon, so order changes if not flippedto flip parameters
            if (extendDirection[1] == 'U' || extendDirection[1] == 'R'){
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

/*! 
  \brief  join two nurbs surfaces prioritising the first surface - as Object3D types

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
Object3D Mapping3D::joinSurfaces(Object3D& srf1, Object3D& srf2, std::string extendDirection){

    cout << "In join surfaces" << endl;

    int deg = srf1.degreeU();// assume the same degree in each direction
    int row_not_col;
    int n_ctrl1;
    int n_ctrl2;
    int n_ctrlCheck2;
    bool newBeforeOld;
    bool flipKnotParam;
    Vector_FLOAT knots1;
    Vector_FLOAT knotsKeep;
    NurbsCurvef crv1;
    NurbsCurvef crv2;
    NurbsCurvef crvOut;
    Vector_HPoint3Df ctrlVec1; // TODO maybe don't make copies here - just use the data directly where possible
    Vector_HPoint3Df ctrlVec2;

    Matrix_HPoint3Df ctrlOut;
    Vector_FLOAT knotsOut;

    Object3D* srfOut;

    cout << "Starting join surfaces" << endl;

    if (srf1.degreeU() != srf2.degreeU() || srf1.degreeV() != srf2.degreeV()){
          cout << "Error: need matching degrees for join direction. Exiting with empty surface" << endl;
          srfOut = new Object3D();
          return *srfOut;// This may break
    }

    // Get whether srf 2 will be using rows (1) or columns (0)
    if (extendDirection[1] == 'L' || extendDirection[1] == 'R'){
        row_not_col = true; // use rows
        n_ctrl2 = srf2.ctrlPnts().cols();       
        n_ctrlCheck2 = srf2.ctrlPnts().rows();  
    }else{
        row_not_col = false; // use columns
        n_ctrl2 = srf2.ctrlPnts().rows();
        n_ctrlCheck2 = srf2.ctrlPnts().cols();  
    }

    cout << "Determine srf 2 extend configuration" << endl;

    // Different cases for rows or columns for srf1
    if (extendDirection[0] == 'L' || extendDirection[0] == 'R'){
        // -----------------------------------------------------------------
        // Row-wise expansion
        // -----------------------------------------------------------------
        // Input checks
        if (n_ctrlCheck2 != srf1.ctrlPnts().rows()){
          cout << "Need to have matched numbers of control points along the join direction\nDoing Knot Insertion" << endl;
          knotInsertionAlongSeam(srf2, extendDirection, srf1.ctrlPnts().rows() - n_ctrlCheck2);
          cout << "Srf1 dimensions are: (" << srf1.ctrlPnts().rows() << ", " << srf1.ctrlPnts().cols() << ")\n";
          cout << "Srf2 dimensions are: (" << srf2.ctrlPnts().rows() << ", " << srf2.ctrlPnts().cols() << ")\n";
        }

        n_ctrl1 = srf1.ctrlPnts().cols();

        cout << "n_ctrl1: " << n_ctrl1 << endl;;

        // Init
        ctrlOut.resize(srf1.ctrlPnts().rows(),n_ctrl1 + n_ctrl2 - deg);
        knotsOut.resize(n_ctrl1 + n_ctrl2 + 1);
              
        knots1 = srf1.knotV();        
        knotsKeep = srf1.knotU();

        cout << "knots1: " << knots1 << endl;
        cout << "knotsKeep: " << knotsKeep << endl;

        // Determine order of combination
        if (extendDirection[0] == 'L'){
            newBeforeOld = true;

            // Determine whether to flip parameters
            if (extendDirection[1] == 'L' || extendDirection[1] == 'D'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }        
        }else{// Extend Right
            newBeforeOld = false;

            // Determine whether to flip parameters
            if (extendDirection[1] == 'R' || extendDirection[1] == 'U'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }

        }
        cout << "Determined order and flipping" << endl;

        // Loop for each row of control points 
        for (int i = 0; i < srf1.ctrlPnts().rows(); i++){
            ctrlVec1 = getMatRow(srf1.ctrlPnts(),i);
            crv1.reset(ctrlVec1,knots1,deg);
            // cout << "reset crv1" << endl;

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

            // cout << "Control vec 2 is: " << ctrlVec2 << endl;
            // cout << "curve 2 (0) is: " << crv2(0) << endl;

            // Join the curves
            crvOut = joinCurves(crv1,crv2,newBeforeOld,flipKnotParam);
            
            // cout << "curve out (0) is: " << crvOut(0) << endl;

            insertMatRow(ctrlOut,crvOut.ctrlPnts(),i);

            knotsOut += crvOut.knot();
        }

        // Average knot vector
        for (int i = 0; i < knotsOut.size(); i++){
            knotsOut[i] /= srf1.ctrlPnts().rows();
        }
        
        // create surface
        cout << "Creating surface" << endl;
        cout << "knots u: " << knotsKeep << " of length: " << knotsKeep.rows() << "\nKnots V: " << knotsOut << " of length: " << knotsOut.rows() << endl;
        cout << "Control points dimensions are: (" << ctrlOut.rows() << ", " << ctrlOut.cols() << ")\n";
        srfOut = new Object3D(deg,deg,knotsKeep,knotsOut,ctrlOut);
        
    }else{
        // -----------------------------------------------------------------
        // Column-wise expansion
        // -----------------------------------------------------------------
        // Input checks
        if (n_ctrlCheck2 != srf1.ctrlPnts().cols()){
          cout << "Need to have matched numbers of control points along the join direction\nDoing Knot Insertion" << endl;
          knotInsertionAlongSeam(srf2, extendDirection, srf1.ctrlPnts().cols() - n_ctrlCheck2);
          cout << "Srf1 dimensions are: (" << srf1.ctrlPnts().rows() << ", " << srf1.ctrlPnts().cols() << ")\n";
          cout << "Srf2 dimensions are: (" << srf2.ctrlPnts().rows() << ", " << srf2.ctrlPnts().cols() << ")\n";
        }

        n_ctrl1 = srf1.ctrlPnts().rows();

        // Init
        ctrlOut.resize(n_ctrl1 + n_ctrl2 - deg, srf1.ctrlPnts().cols());
        knotsOut.resize(n_ctrl1 + n_ctrl2 + 1);

        knots1 = srf1.knotU(); 
        knotsKeep = srf1.knotV();

        cout << "knots1: " << knots1 << endl;
        cout << "knotsKeep: " << knotsKeep << endl;       

        // Determine order of combination
        if (extendDirection[0] == 'D'){
            newBeforeOld = true;

            // Determine whether to flip parameters
            if (extendDirection[1] == 'D' || extendDirection[1] == 'L'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }        
        }else{// Extend Up
            newBeforeOld = false;

            // Determine whether // take columns in same order: assume a rotati// take columns in same order: assume a rotation, so order changes if not flippedon, so order changes if not flippedto flip parameters
            if (extendDirection[1] == 'U' || extendDirection[1] == 'R'){
                flipKnotParam = false;
            }else{
                flipKnotParam = true;
            }

        }

        // Loop for each column of control points 
        for (int j = 0; j < srf1.ctrlPnts().cols(); j++){
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
                    // take points in the columns in different order - from the end first 
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
            knotsOut[i] /= srf1.ctrlPnts().cols();
        }
        
        // create surface
        cout << "Creating surface" << endl;
        cout << "knots V: " << knotsKeep << " of length: " << knotsKeep.rows() << "\nKnots U: " << knotsOut << " of length: " << knotsOut.rows() << endl;
        cout << "Control points dimensions are: (" << ctrlOut.rows() << ", " << ctrlOut.cols() << ")\n";
        srfOut = new Object3D(deg,deg,knotsOut,knotsKeep,ctrlOut);

    }
    

    return *srfOut;

}

/*! 
  \brief  Compute the number of control points required for an extension of a surface

  \param extendDirection string to indicate how the surfaces are to be joined
  Options are L, R, U, D. Is a string of two values, one for srf1 and one for srf2
  \param data             the data on the new extension
  \param srf              The surface being extended

  \author Benjamin Morrell
  \date 28 March 2018
*/
std::vector<int> Mapping3D::computeNumberOfControlPoints(std::string extendDirection, Matrix_Point3Df& data, Object3D& srf){
  // TDBM - later - just input the object ID and take the object from the map
  
  // Standard number of control points to use
  int nCtrlStandard = std::min(srf.ctrlPnts().rows(),srf.ctrlPnts().cols());

  cout << "Srf control points: " << srf.ctrlPnts().rows() << " rows, " << srf.ctrlPnts().cols() << " cols.\n";

  cout << "Frac 1: " << (float)data.cols()/(float)data.rows() << endl;
  cout << "Frac 2: " << (float)data.rows()/(float)data.cols() << endl;

  // Init output
  std::vector<int> nCtrlNew(2);


  // Match the number of control points in the direction it will be joined
  if (extendDirection[1] == 'L' || extendDirection[1] == 'R'){
    // It will join along the t direction, need to match the s number
    if (extendDirection[0] == 'L' || extendDirection[0] == 'R' ){
      // Match the s direction to the s direction
      nCtrlNew[0] = srf.ctrlPnts().rows();
      //Set t direction to be a fraction of the standard number (assume less cols than rows)
      nCtrlNew[1] = (int)(round((float)data.cols()/(float)data.rows()*nCtrlStandard));
      // TDBM check this gives the correct output
    }else{
      //Match the s direction to the t direction
      nCtrlNew[0] = srf.ctrlPnts().cols();
      //Set t direction to be a fraction of the standard number (assume less cols than rows)
      nCtrlNew[1] = (int)(round((float)data.cols()/(float)data.rows()*nCtrlStandard));
    }
  }else if (extendDirection[1] == 'D' || extendDirection[1] == 'U'){
    // It will join along the s direction, need to match the t number
    if (extendDirection[0] == 'L' || extendDirection[0] == 'R' ){
      // Match the s direction to the s direction
      nCtrlNew[1] = srf.ctrlPnts().rows();
      //Set t direction to be a fraction of the standard number (assume less cols than rows)
      nCtrlNew[0] = (int)(round((float)data.rows()/(float)data.cols()*nCtrlStandard));
      // TDBM check this gives the correct output
    }else{
      //Match the s direction to the t direction
      nCtrlNew[1] = srf.ctrlPnts().cols();
      //Set t direction to be a fraction of the standard number (assume less cols than rows)
      nCtrlNew[0] = (int)(round((float)data.rows()/(float)data.cols()*nCtrlStandard));
    }
  }else{
    nCtrlNew[0] = 0;
    nCtrlNew[1] = 0;
  }

  

  // Increase number of knots in join direction to be above 2, if criteria is met
  if (nCtrlNew[0] == 2 || nCtrlNew[0] == 3){
    // if ((float)data.rows()/(float)data.cols()*srf.ctrlPnts().rows() > 1.8){
      nCtrlNew[0] = 4;
    // }
  }else if (nCtrlNew[1] == 2 || nCtrlNew[1] == 3){
    // if ((float)data.cols()/(float)data.rows()*srf.ctrlPnts().cols() > 1.8){
      nCtrlNew[1] = 4;
    // }
  }

  // Make sure there are not more control points than data
  if (nCtrlNew[0] > data.rows()){
    cout << "\n\tMore control points: " << nCtrlNew[0] << " than data: " << data.rows() << ". Reducing control points\n" << endl;
    nCtrlNew[0] = data.rows();
  }
  if (nCtrlNew[1] > data.cols()){
    cout << "\n\tMore control points: " << nCtrlNew[1] << " than data: " << data.cols() << ". Reducing control points\n" << endl;
    nCtrlNew[1] = data.cols();
  }

  return nCtrlNew;

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
void Mapping3D::meshFromScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloudOut, pcl::PointCloud<pcl::PointNormal>::Ptr cloudIn){
  // Initialise
  Eigen::Array<bool, Eigen::Dynamic, 1> rowFlags(cloudIn->height, 1); rowFlags.setOnes(cloudIn->height, 1);
  Eigen::Array<bool, 1, Eigen::Dynamic> colFlags(1, cloudIn->width); colFlags.setOnes(1, cloudIn->width);
  // Initialise Nan array
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> nanArray(cloudIn->height,cloudIn->width);
  nanArray.setZero(cloudIn->height,cloudIn->width);
  bRejectScan = false;

  std::cerr << "PointCloud before filtering: W: " << cloudIn->width << "\tH: " << cloudIn->height << "\tdata points." << std::endl;

  // Fill the Nan Array
  getNanMatrixFromPointCloud(nanArray, cloudIn);

  cout << "Number of Nans: " << nanArray.count() << endl;

  if ((float)nanArray.count()/(float)(cloudIn->width*cloudIn->height) > 0.98){
    cout << "Too many Nans in the scan. Rejecting" << endl;
    bRejectScan = true;
    return;
  } 

  bool exitFlag = false;

  bool nansPresent = true;

  // Loop to remove Nans
  while (!exitFlag){
    nansPresent = removeRowsNan(nanArray, rowFlags);

    if (nansPresent){
      nansPresent = removeColsNan(nanArray, colFlags);
      if (nansPresent){
        // Check dimensions 
        // cout << "Nans removed: Row count: " << rowFlags.count() << "\nCol Count: " << colFlags.count() << endl;
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

  cout << "Number of nans left: " << nanArray.count() << endl;


  if ((float)nanArray.count()/(float)(numRowsDesired*numColsDesired) > 0.8){
    cout << "Too many Nans in the scan after downsampling. Rejecting" << endl;
    bRejectScan = true;
    return;
  } 

  cout << "Number of rows: " << rowFlags.count() << endl;
  cout << "Number of cols: " << colFlags.count() << endl;

  int numRows = this->numRowsDesired;
  int numCols = this->numColsDesired;

  // Rest cloud size given number of output rows and cols
  if (rowFlags.count() < numRowsDesired){
    pcl::common::deleteRows(*cloudOut, *cloudOut, (int)std::ceil((float)(numRowsDesired - rowFlags.count())/2.0));
    // pcl::common::deleteRows(*cloudOut, *cloudOut, numRowsDesired - rowFlags.count());
    cout << "Remove rows from output. Size is now: Rows: " << cloudOut->height << ", Cols: " << cloudOut->width << endl;
    // Note that delete Rows removes twice the input amount -= so take the new values from the output cloud
    numRows = cloudOut->height;
  }

  if (colFlags.count() < numColsDesired){
    pcl::common::deleteCols(*cloudOut, *cloudOut, (int)std::ceil((float)(numColsDesired - colFlags.count())/2.0));
    // pcl::common::deleteCols(*cloudOut, *cloudOut, numColsDesired - colFlags.count());
    cout << "Remove cols from output. Size is now: Rows: " << cloudOut->height << ", Cols: " << cloudOut->width << endl;
    // Note that delete Cols removes twice the input amount -= so take the new values from the output cloud
    numCols = cloudOut->width;
  }

  // Downsample 
  downsampleRow(rowFlags, numRows);
  downsampleCol(colFlags, numCols);

  cout << "After downsample:" << endl;
  cout << "Number of rows: " << rowFlags.count() << endl;
  cout << "Number of cols: " << colFlags.count() << endl;
  cout << "Cloud size is: " << cloudOut->height << ", Cols: " << cloudOut->width << endl;

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
          cloudOut->at(jj,ii).x = cloudIn->at(j,i).x;
          cloudOut->at(jj,ii).y = cloudIn->at(j,i).y;
          cloudOut->at(jj,ii).z = cloudIn->at(j,i).z;

          
          // Store indices if the value is nan
          if (!pcl::isFinite(cloudIn->at(j,i))){
            nanIndices(0,ijk) = ii;
            nanIndices(1,ijk) = jj;
            ijk++; 
          }
          // else if (std::abs(cloudOut->at(jj,ii).z) < 3.0 || std::abs(cloudOut->at(jj,ii).z) > 10.0){
          //   // Filter out large z variances THESE ARE HARD CODED FOR THE MOMENT!!!
          //   cout << "Removing point because out of z threshold. z was: " << cloudOut->at(jj,ii).z << endl;
            
          //   cloudOut->at(jj,ii).x = 1.0/0.0;
          //   cloudOut->at(jj,ii).y = 1.0/0.0;
          //   cloudOut->at(jj,ii).z = 1.0/0.0;

          //   nanIndices(0,ijk) = ii;
          //   nanIndices(1,ijk) = jj;
          //   ijk++; 
          // }

          // Increment the new column index
          jj++;
        }
      }
      // Increment the new row index
      ii++;
    }
  }

  cout << "Finished copying data across" << endl;

  cout << "Number of nans left: " << nanArray.count() << endl;

  // pcl::PCDWriter writer;
  // writer.write<pcl::PointNormal> ("cloud_pre_nan_average.pcd", *cloudOut, false);

  // cout << "NanIndices are: " << nanIndices << endl;

  // Average Nans
  bool noNans = false;
  while (!noNans){
    noNans = averageOutNans(cloudOut, nanIndices);
    cout << "Iteration for averaging Nans" << endl;
  }

  // writer.write<pcl::PointNormal> ("cloud_post_nan_average.pcd", *cloudOut, false);

}

/*! 
  \brief  Searches through an organised pointcloud and makes a boolean array with 1s where the point cloud is nan

  \param nanArray   the output boolean array
  \param cloud      a pointer to the point cloud to process
  
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::getNanMatrixFromPointCloud(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray, pcl::PointCloud<pcl::PointNormal>::Ptr cloud){

  float zLimLow = -99999.9;
  float zLimHigh = 99999.9;


  if (bFilterZ){
    // Compute mean and covariance
    std::vector<int> indices;
    int nPoints = nPointsZLim; // TODO - make this a parameter
    int nData = cloud->width*cloud->height;
    for (int i = 0; i < nData; i = i + nData/nPoints){indices.push_back(i);}

    Eigen::Matrix<float, 3, 3> covarianceMatrix;
    Eigen::Matrix<float, 4, 1> centroid;

    pcl::computeMeanAndCovarianceMatrix(*cloud, indices, covarianceMatrix, centroid);

    cout << "Covariance matrix is: " << covarianceMatrix << endl;
    cout << "Centroid is: " << centroid << endl;
    // TODO reuse these for computing metrics?

    float threeSigZ = std::sqrt(covarianceMatrix(2,2))*3.0;

    zLimLow = centroid(2) - threeSigZ;
    zLimHigh = centroid(2) + threeSigZ*zThreshMultiplier;
  }

  cout << "Zlow is: " << zLimLow << ", Zhigh is: " << zLimHigh << endl;

   // Loop through Point cloud 
  for (int i = 0; i < cloud->height; i++){
    for (int j = 0; j < cloud->width; j++){
      // If nan value
      if (!pcl::isFinite(cloud->at(j,i))){ // cloud->at(col,row)
        nanArray(i,j) = true;
      }else if (cloud->at(j,i).z > zLimHigh){
        nanArray(i,j) = true;
        // cout << "Capping Z high at value: " << cloud->at(j,i).z << endl;
      }else if (cloud->at(j,i).z < zLimLow){
        nanArray(i,j) = true;
        // cout << "Capping Z low at value: " << cloud->at(j,i).z << endl;
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

  // cout << "Max NaN count: " << maxNan << "\tmax allowed is: " << this->maxNanAllowed << endl;

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

  // cout << "Max NaN count: " << maxNan << "\tmax allowed is: " << this->maxNanAllowed << endl;

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
void Mapping3D::downsampleRow(Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags, int numRows){

  // Flag for default input
  if (numRows < 0 ){
    numRows = this->numRowsDesired;
  }

  // Exit if already few enough rows or columnts
  if (rowFlags.count() <= numRows){
    // Nothing to do
    cout << "Number of rows less than or equal to desired. Not downsampling" << endl;
    return;
  }

  // Create linspaced array 
  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArrayf.setLinSpaced(numRows, 0, rowFlags.count()-1);
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
void Mapping3D::downsampleCol(Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags, int numCols){

  // Flag for default input
  if (numCols < 0 ){
    numCols = this->numColsDesired;
  }

  // Exit if already few enough rows or columnts
  if (colFlags.count() <= numCols){
    // Nothing to do
    cout << "Number of cols less than or equal to desired. Not downsampling" << endl;
    return;
  }

  // Create linspaced array 
  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArrayf.setLinSpaced(numCols, 0, colFlags.count()-1);
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
bool Mapping3D::averageOutNans(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic>& nanIndices){
  // cout << "Number of cols: " << nanIndices.cols() << endl;

  bool noNans = true;

  for (int i = 0; i < nanIndices.cols(); i++){
    if (nanIndices(0,i) == -1){
      // Flag meaning that those terms are not valid (not used)
      break;
    }else if (nanIndices(0,i) == -2){
      // Flag meaning that those terms are not valid (not used), but are not at the end
      continue;
    }
    // cout << "Before: " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
    regionAverageSimple(cloud,nanIndices(0,i),nanIndices(1,i));
    // cout << "After:  " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
    if (!pcl::isFinite(cloud->at(nanIndices(1,i),nanIndices(0,i)) )){
      // Still nan, keep as flag
      noNans = false;
    }else{
      // Indicate they are not active
      nanIndices(0,i) = nanIndices(1,i) = -2;
    }
  }

  cout << "All nans gone?: " << noNans << endl;

  return noNans;

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
void Mapping3D::regionAverage(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, int i, int j){

  pcl::PointXYZ average(0.0,0.0,0.0);

  int add_count = 0; // to track how many values are averaged. 

  // Back in i
  if (i > 0){
    if (pcl::isFinite(cloud->at(j,i-1))){
      if (i > 1){
        if (pcl::isFinite(cloud->at(j,i-2))){
          // Step with same delta as neighbours
          // average.getArray3fMap() += 2*cloud->at(j,i-1).getArray3fMap() - cloud->at(j,i-2).getArray3fMap();
          average.x += 2*cloud->at(j,i-1).x - cloud->at(j,i-2).x;
          average.y += 2*cloud->at(j,i-1).y - cloud->at(j,i-2).y;
          average.z += 2*cloud->at(j,i-1).z - cloud->at(j,i-2).z;
        }else{
          // average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
          average.x += cloud->at(j,i-1).x;
          average.y += cloud->at(j,i-1).y;
          average.z += cloud->at(j,i-1).z;
        }

      }else{
        // average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
        average.x += cloud->at(j,i-1).x;
        average.y += cloud->at(j,i-1).y;
        average.z += cloud->at(j,i-1).z;        
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
          // average.getArray3fMap() += 2*cloud->at(j,i+1).getArray3fMap() - cloud->at(j,i+2).getArray3fMap();   
          average.x += 2*cloud->at(j,i+1).x - cloud->at(j,i+2).x;
          average.y += 2*cloud->at(j,i+1).y - cloud->at(j,i+2).y;
          average.z += 2*cloud->at(j,i+1).z - cloud->at(j,i+2).z; 
        }else{
          // average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
          average.x += cloud->at(j,i+1).x;
          average.y += cloud->at(j,i+1).y;
          average.z += cloud->at(j,i+1).z;
        }
      }else{
        // average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
        average.x += cloud->at(j,i+1).x;
        average.y += cloud->at(j,i+1).y;
        average.z += cloud->at(j,i+1).z;
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
          // average.getArray3fMap() += 2*cloud->at(j-1,i).getArray3fMap() - cloud->at(j-2,i).getArray3fMap();
          average.x += 2*cloud->at(j-1,i).x - cloud->at(j-2,i).x;
          average.y += 2*cloud->at(j-1,i).y - cloud->at(j-2,i).y;
          average.z += 2*cloud->at(j-1,i).z - cloud->at(j-2,i).z;
        }else{
          // average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
          average.x += cloud->at(j-1,i).x;
          average.y += cloud->at(j-1,i).y;
          average.z += cloud->at(j-1,i).z;
        }
      }else{
        // average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
        average.x += cloud->at(j-1,i).x;
        average.y += cloud->at(j-1,i).y;
        average.z += cloud->at(j-1,i).z;
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
          // average.getArray3fMap() += 2*cloud->at(j+1,i).getArray3fMap() - cloud->at(j+2,i).getArray3fMap();
          average.x += 2*cloud->at(j+1,i).x - cloud->at(j+2,i).x;
          average.y += 2*cloud->at(j+1,i).y - cloud->at(j+2,i).y;
          average.z += 2*cloud->at(j+1,i).z - cloud->at(j+2,i).z;
        }else{
          // average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
          average.x += cloud->at(j+1,i).x;
          average.y += cloud->at(j+1,i).y;
          average.z += cloud->at(j+1,i).z;  
        }
      }else{
        // average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
        average.x += cloud->at(j+1,i).x;
        average.y += cloud->at(j+1,i).y;
        average.z += cloud->at(j+1,i).z;        
      }
      add_count ++;
    }    
  }

  // Divide by count to get the average
  average.getArray3fMap() /= (float)add_count;

  //update cloud
  cloud->at(j,i).x = average.x;
  cloud->at(j,i).y = average.y;
  cloud->at(j,i).z = average.z;

}

/*! 
  \brief  Replaces a given point by averaginf from neighbouring points

  \param cloud  pointer to cloud to process
  \param i      row index of point to change
  \param j      col index of point to change
  
  Will search +_ in i and j, then try to take the average of them
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::regionAverageSimple(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, int i, int j){

  pcl::PointXYZ average(0.0,0.0,0.0);

  int add_count = 0; // to track how many values are averaged. 

  // Back in i
  if (i > 0){
    if (pcl::isFinite(cloud->at(j,i-1))){
      // average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
      average.x += cloud->at(j,i-1).x;
      average.y += cloud->at(j,i-1).y;
      average.z += cloud->at(j,i-1).z;        
      add_count ++;
    }    
  }

  // Forward in i
  if (i < cloud->height-1){
    if (pcl::isFinite(cloud->at(j,i+1))){
        // average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
        average.x += cloud->at(j,i+1).x;
        average.y += cloud->at(j,i+1).y;
        average.z += cloud->at(j,i+1).z;
      add_count ++;
    }
  }

  // Back in j
  if (j > 0){
    if (pcl::isFinite(cloud->at(j-1,i))){
      // average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
      average.x += cloud->at(j-1,i).x;
      average.y += cloud->at(j-1,i).y;
      average.z += cloud->at(j-1,i).z;
      add_count ++;
    }    
  }

  // Forward in j
  if (j < cloud->width-1){
    if (pcl::isFinite(cloud->at(j+1,i))){
      // average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
      average.x += cloud->at(j+1,i).x;
      average.y += cloud->at(j+1,i).y;
      average.z += cloud->at(j+1,i).z;
      add_count ++;
    }
  }

  // Divide by count to get the average
  average.getArray3fMap() /= (float)add_count;

  //update cloud
  cloud->at(j,i).x = average.x;
  cloud->at(j,i).y = average.y;
  cloud->at(j,i).z = average.z;

}

//-------------------------------------------------------------------
// ------------- UPDATE OBJECT  -----------------------------
//-------------------------------------------------------------------
/*! 
  \brief  Adds an object to the map

  \param cloud          the point cloud for the new object, in the global frame
  \param searchMetrics  the metrics describing the new object
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
void Mapping3D::addObject(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, std::vector<float> searchMetrics){

  // cout << "In addObject, cloud at (0,0) is: " << cloud->at(0,0) << endl;
  // Convert data to nurbs data
  Matrix_Point3Df mesh = nurbsDataFromPointCloud(cloud);

  // cout << "Mesh from cloud point (0,0) is: " << mesh(0,0) << endl;
  // cout << "Mesh from cloud point (5,5) is: " << mesh(5,5) << endl;
  cout << "nCtrlDefault is: " << nCtrlDefault[0] << ", " << nCtrlDefault[1] << endl;

  // Compute number of control points (TDBM  - imrpove this)
  int nCtrl1 = std::min((int)cloud->height/2,nCtrlDefault[0]);
  int nCtrl2 = std::min((int)cloud->width/2,nCtrlDefault[1]);
  cout << "Number of control points computed is: (" << nCtrl1 << ", " << nCtrl2 << ")\n";
  Object3D obj(mesh, order[0], order[1], nCtrl1, nCtrl2);

  // cout << "object point (0,0) " << obj(0.0,0.0) << endl;
  // cout << "object point (0.4,0.4) " << obj(0.4,0.4) << endl;

  // Add NURBS Object (using default settings)
  objectMap.push_back(obj);

  // TDBM - consider whether to do this or not... NOTE that this changes from data association
  searchMetrics[0] = objectMap[objectMap.size()-1].getCentre().x();
  searchMetrics[1] = objectMap[objectMap.size()-1].getCentre().y();
  searchMetrics[2] = objectMap[objectMap.size()-1].getCentre().z();
  searchMetrics[3] = objectMap[objectMap.size()-1].getObjSize();

  // check search metrics size

  // Add search metrics
  objectMetrics.push_back(searchMetrics);

  cout << "New object added with ID: " << objectMap.size()-1 << " and with metrics: ";
  for (int i = 0; i < numberOfMetrics; i++){cout << searchMetrics[i] << ", ";}cout << endl;

}

/*! 
  \brief  Adds an object from file 

  \param objID  id of the object to update
  \param obj    Object3D parameters to replace those in the map
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
void Mapping3D::addObjectFromFile(const char * filename){

  // Checks if filename exists - TODO

  // Load from filename
  Object3D obj;

  obj.readObject3D(filename);

  // Add object in objectMap
  objectMap.push_back(obj);

  // Compute new metrics
  std::vector<float> searchMetrics = computeSearchMetrics(obj);

  // Add metrics
  objectMetrics.push_back(searchMetrics);

  cout << "New object added with ID: " << objectMap.size()-1 << " and with metrics: ";
  for (int i = 0; i < numberOfMetrics; i++){cout << searchMetrics[i] << ", ";}cout << endl;

}


/*! 
  \brief  Updates an object in the map

  \param objID  id of the object to update
  \param obj    Object3D parameters to replace those in the map
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
void Mapping3D::updateObjectInMap(int objID, Object3D& obj){

  if (objID > objectMap.size()-1){
    cout << "Error: input ID larger than object map. Not updating\n\n" << endl;
    return;
  }

  // replace object in objectMap
  objectMap[objID] = obj;

  // Compute new metrics
  std::vector<float> searchMetrics = computeSearchMetrics(obj);

  // Update metrics
  objectMetrics[objID] = searchMetrics;

  cout << "Object " << objID << " updated with metrics: ";
  for (int i = 0; i < numberOfMetrics; i++){cout << searchMetrics[i] << ", ";}cout << endl;

}


/*! 
  \brief  Updates an object in the map with a new observation

  \param objID    id of the object to update
  \param obsObjPC point cloud of the new observataion 
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
void Mapping3D::updateObject(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC){
  
  // Get point cloud from NURBS object
  pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC(new pcl::PointCloud<pcl::PointNormal>(mtSurf, msSurf, pcl::PointNormal()));
  pointCloudFromObject3D(objID, msSurf, mtSurf, mapObjPC);

  // cout << "Point cloud from object at (0,0) is: " << mapObjPC->at(0,0) << endl;
  // cout << "Point cloud from observation at (0,0) is: " << obsObjPC->at(0,0) << endl;
  // // Save point clouds to Matlab
  // // Write the downsampled version to disk
  // pcl::PCDWriter writer;
  // writer.write<pcl::PointNormal> ("observation.pcd", *obsObjPC, false);
  // writer.write<pcl::PointNormal> ("map_data.pcd", *mapObjPC, false);

  
  // Split new surface observation
  SplitSurface ss;
  ss.newRowColBuffer = newRowColBuffer; // number of allowed overlap points in a row/column

  ss.splitNewSurfaceObservation(mapObjPC, obsObjPC);
  cout << "Extend Direction is: " << ss.extendDirection << endl;
  

  if (ss.extendDirection[0] == 'N' || ss.extendDirection[1] == 'N'){
    cout << "\nNo extension of surface needed (not enough new data)\n" << endl;
    return;
  }
  

  if (useNonRectData){
    // Get new data indices 
    ss.getNewDataIndices();
    cout << "New Row indices are:\n" << ss.newRowIndices << endl;
    cout << "New Col indices are:\n" << ss.newColIndices << endl;
  }

  // Get matrix from new data
  Matrix_Point3Df mesh;
  if (useNonRectData){
    mesh = nurbsDataFromPointCloud(obsObjPC, ss.newRowIndices, ss.newColIndices);
  }else{
    mesh = nurbsDataFromPointCloud(obsObjPC, ss.newDataIndices);
  }
  cout << "Mesh2 size is (" << mesh.rows() << ", " << mesh.cols() << ")\n";

  // Get the number of control points
  std::vector<int> nCtrlNew = computeNumberOfControlPoints(ss.extendDirection, mesh, objectMap[objID]);
  cout << "new control points: " << nCtrlNew[0] << ", " << nCtrlNew[1] << endl;

  if (nCtrlNew[0] <= order[0] || nCtrlNew[1] <= order[1]){
    cout << "\nNot enough control points, no extension\n" << endl;
    return;
  }

  // Create Object for new surface? (nknots = deg + nctrl + 1)
  Object3D obj(mesh, order[0], order[1], nCtrlNew[0], nCtrlNew[1]);

  // if (useNonRectData){
  //   // CHECK OBJECT:
  obj.writeVRML("newSurf.wrl",Color(255,100,255),50,80);

  // Get point cloud from NURBS object
  pcl::PointCloud<pcl::PointNormal>::Ptr objPC(new pcl::PointCloud<pcl::PointNormal>(10, 10, pcl::PointNormal()));
  obj.getSurfacePointCloud(objPC, 10, 10);
  pcl::PCDWriter writer;
  writer.write<pcl::PointNormal> ("newSurf.pcd", *objPC, false);

  //   // Test generating data
  //   int nPoints = 25;
    
  //   // Get the data points
  //   Matrix_Point3Df data = obj.getSurfacePoints(nPoints, nPoints);

  //   bool bIsnoNan = true;
  //   for (int i = 0; i < nPoints; i++){
  //     for (int j = 0; j < nPoints; j++){
  //       if (!std::isfinite(data(i,j).x())){
  //           bIsnoNan = false;
  //           cout << "\n\n\t\tNANS IN New object\n\n" << endl;
  //       }
  //     }
  //   }
  // }
 
  
  // cout << "object of new mesh has data at (0,0): " << obj.pointAt(0.0,0.0) << endl;  
  // cout << "object of new mesh has data at (1.0,1.0): " << obj.pointAt(1.0,1.0) << endl;

  // Knot insertion
  if (knotInsertionFlag){
    cout << "inserting knots" << endl;
    knotInsertionPreMerge(obj, ss.extendDirection);
  }
  
  // Get point cloud from NURBS object
  // obj.getSurfacePointCloud(objPC, 10, 10);
  // writer.write<pcl::PointNormal> ("newSurf_knot_insert.pcd", *objPC, false);


  // Join Surfaces
  cout << "joining surfaces" << endl;
  Object3D objJoin = joinSurfaces(objectMap[objID],obj,ss.extendDirection);

  cout << "Point (0.3,0.3) on obj3: " << objJoin.pointAt(0.3,0.3) << endl;

  // Update the map
  updateObjectInMap(objID, objJoin);

  objectMap[objID].writeVRML("updatedSurf.wrl",Color(255,100,255),50,80);

  // Get point cloud from NURBS object
  writeObjectPCDFile("updatedSurf.pcd", objID, 10, 10);
  
}

/*! 
  \brief  inserts knots into a surface in preparation for merging with another surface

  \param      obj           the Object3D to insert knots for
  \param[in]  exendDirection The direction surfaces are joining
  
  \author Benjamin Morrell
  \date 5 April 2018
*/
void Mapping3D::knotInsertionPreMerge(Object3D& obj, std::string extendDirection){

  // From global settings
  // float deltaKnotInsert = 1e-4;
  // int numInsert = 3;

  Vector_FLOAT knotsInsert(numInsert);

  int nCtrlS = obj.ctrlPnts().rows();
  int nCtrlT = obj.ctrlPnts().cols();

  float startVal;
  float endVal;
  bool insertUFlag;

  switch (extendDirection[1]){
    case 'L' :
      // Add knots in the t direction near the right end (end of parameters)
      startVal = obj.knotV(nCtrlT-1) + deltaKnotInsert;
      endVal = obj.knotV(nCtrlT) - deltaKnotInsert;
      insertUFlag = false;
      break;
    case 'R' : 
      // Add knots in the t direction near the left end (start of parameters)
      startVal = obj.knotV(obj.degreeV()) + deltaKnotInsert;
      endVal = obj.knotV(obj.degreeV()+ 1) - deltaKnotInsert;
      insertUFlag = false;
      break;
    case 'D' :
      // Add knots in the s direction near the upper end (end of parameters)
      startVal = obj.knotU(nCtrlS-1) + deltaKnotInsert;
      endVal = obj.knotU(nCtrlS) - deltaKnotInsert;
      insertUFlag = true;
      break;
    case 'U' : 
      // Add knots in the s direction near the lower end (start of parameters)
      startVal = obj.knotU(obj.degreeU()) + deltaKnotInsert;
      endVal = obj.knotU(obj.degreeU() + 1) - deltaKnotInsert;
      insertUFlag = true;
      break;
  }

  // Fill knot vector
  float step = (endVal - startVal)/((float)numInsert-1);

  for (int i = 0; i < numInsert; i++){
    knotsInsert[i] = startVal + step * i;
  }

  

  // Do knot insertion
  if (insertUFlag){
    cout << "knots pre-insertion are: " << obj.knotU() << endl;
    cout << "Inserting knots in U: " << knotsInsert << endl;
    obj.refineKnotU(knotsInsert);
  }else{
    cout << "knots pre-insertion are: " << obj.knotV() << endl;
    cout << "Inserting knots in V: " << knotsInsert << endl;
    obj.refineKnotV(knotsInsert);
  }
}

/*! 
  \brief  inserts knots into a surface in preparation for merging with another surface - along the seam

  \param      obj           the Object3D to insert knots for
  \param[in]  exendDirection The direction surfaces are joining
  \param[in]  nInsert         The number of knots to insert
  
  \author Benjamin Morrell
  \date 5 April 2018
*/
void Mapping3D::knotInsertionAlongSeam(Object3D& obj, std::string extendDirection, int nInsert){

  // From global settings
  // float deltaKnotInsert = 1e-4;
  // int numInsert = 3;

  Vector_FLOAT knotsInsert(nInsert);

  int nCtrlS = obj.ctrlPnts().rows();
  int nCtrlT = obj.ctrlPnts().cols();

  float startVal = 0.0 + deltaKnotInsert;
  float endVal = 1.0 - deltaKnotInsert;
  bool insertUFlag;

  switch (extendDirection[1]){
    case 'L' :
      // Add knots in the s direction 
      insertUFlag = true;
      break;
    case 'R' : 
      // Add knots in the s direction 
      insertUFlag = true;
      break;
    case 'D' :
      // Add knots in the t direction 
      insertUFlag = false;
      break;
    case 'U' : 
      // Add knots in the t direction 
      insertUFlag = false;
      break;
  }

  // Fill knot vector
  float step = (endVal - startVal)/((float)nInsert-1);

  for (int i = 0; i < nInsert; i++){
    knotsInsert[i] = startVal + step * i;
  }

  

  // Do knot insertion
  if (insertUFlag){
    cout << "AlongSeam knots pre-insertion are: " << obj.knotU() << endl;
    cout << "AlongSeam Inserting knots in U: " << knotsInsert << endl;
    obj.refineKnotU(knotsInsert);
  }else{
    cout << "AlongSeam knots pre-insertion are: " << obj.knotV() << endl;
    cout << "AlongSeam Inserting knots in V: " << knotsInsert << endl;
    obj.refineKnotV(knotsInsert);
  }
}

//-------------------------------------------------------------------
// ------------- DATA ASSOCIATION  -----------------------------
//-------------------------------------------------------------------

/*! 
  \brief  Computes search metrics from an Object3D object

  \param[in]  obj           the Object3D to compute the metrics for
  \param[out] searchMetrics the list of metrics for data association searches
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
std::vector<float> Mapping3D::computeSearchMetrics(Object3D& obj){
  
  std::vector<float> searchMetrics(numberOfMetrics);

  // Centre metrics - take from object
  searchMetrics[0] = obj.getCentre().x();
  searchMetrics[1] = obj.getCentre().y();
  searchMetrics[2] = obj.getCentre().z();

  // Size
  searchMetrics[3] = obj.getObjSize();

  // Color
  searchMetrics[4] = obj.getColor().x();
  searchMetrics[5] = obj.getColor().y();
  searchMetrics[6] = obj.getColor().z();

  return searchMetrics;
}

/*! 
  \brief  Computes search metrics from a pcl point cloud

  \param[in]  cloud         the pcl point cloud to compute the metrics for
  \param[out] searchMetrics the list of metrics for data association searches
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
std::vector<float> Mapping3D::computeSearchMetrics(pcl::PointCloud<pcl::PointNormal>::Ptr cloud){
  
  // cout << "In ComputeSearchMetrics\nPC(0,0) is: " << cloud->at(0,0) << endl;

  // init
  std::vector<float> searchMetrics(numberOfMetrics);

  // Get Centroid
  Eigen::Vector4f centroid; // Last component is 1 to allow use of a 4x4 transformation matrix
  pcl::compute3DCentroid(*cloud,centroid); 

  // Compute size (maximum x, y, z dimension) TDBM make this a better measure
  pcl::PointNormal minP;
  pcl::PointNormal maxP;
  pcl::getMinMax3D(*cloud, minP, maxP);

  // cout << "after centroid and minmax, PC(0,0) is: " << cloud->at(0,0) << endl;

  float objSize = std::max(maxP.x - minP.x,std::max(maxP.y - minP.y,maxP.z - minP.z));

  // Centre metrics
  searchMetrics[0] = centroid.x();
  searchMetrics[1] = centroid.y();
  searchMetrics[2] = centroid.z();

  // Size
  searchMetrics[3] = objSize;

  // Color
  searchMetrics[4] = 0.0;
  searchMetrics[5] = 0.0;
  searchMetrics[6] = 0.0;

  return searchMetrics;
}

/*! 
  \brief  Use input search metrics to find a matching object in the objectMap

  \param searchMetrics vector of metrics to be searched for

  uses: objectMetrics
        numberOfMetrics
        searchThresh
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
int Mapping3D::dataAssociation(std::vector<float> searchMetrics){

  // Input
  // Eigen::Array<float, 1, Eigen::Dynamic> searchMetrics(1,7);

  // From object map [x, y, z, size, r, g, b]. 
  // std::vector<float> searchThresh;

  int objID = -1;
  float dist;

  if (objectMap.size() == 0){
    // No objects in the map
    return objID;
  }

  // Distances array
  Eigen::Array<float, 1, Eigen::Dynamic> distances(1,objectMap.size());
  distances.setZero(objectMap.size());

  Eigen::Array<bool, 1, Eigen::Dynamic> activeObjects(1,objectMap.size());
  activeObjects.setOnes(objectMap.size()); 

  // Search through each metric
  for (int j = 0; j < numberOfMetrics; j++){
    if (activeObjects.count() < 1){
      // No match 
      cout << "No Object Match - new object." << endl;
      objID = -1;
      return objID;
    }
    if (searchThresh[j] > 0.0f){
      // If search metric is active

      // Get distance for each object
      for (int i = 0; i < objectMetrics.size(); i ++){
        if (activeObjects[i]){// If the object is still active
          // Distance for that metric
          dist = pow(objectMetrics[i][j]-searchMetrics[j],2.0f);

          // Check against threshold
          if (dist < searchThresh[j]){
            // Add distances to get cumulative distance
            distances[i] += dist;
          }else{
            // Fail criteria form search
            activeObjects[i] = false;
            cout << "Object " << i << " fails search criteria" << endl;
            // Set distance to high to assist search
            // distances[i] = 999999.9;
          }
        }
      }
    }
  }

  // If there are multiple matches - get the min distance
  int i;
  float minDist = 99999.9;
  if (activeObjects.count() == 1){
    // Get the object ID
    cout << "Only one object in search criteria" << endl;
    activeObjects.maxCoeff(&objID);
  }else if (activeObjects.count() < 1){
    cout << "No objects in search criteria" << endl;
    objID = -1;
  }else{
    cout << "Multiple objects in search criteria" << endl;
    while (activeObjects.count() > 0){
      // Get index of the first remaining match
      activeObjects.maxCoeff(&i);
      
      // Get the distance 
      dist = distances[i];

      if (dist < minDist){
        // Select objID to get min
        objID = i;
        // update min dist
        minDist = dist;
        cout << "min dist is " << minDist << "with id " << i << endl;
      }
      // Update activeObjects to get next search
      activeObjects[i] = false;
    }
  }
  cout << "Object distances are: " << distances << endl;
  cout << "Object match is: " << objID << endl;

  return objID;

  // There may be a more efficient data type to use than an Eigen Array - resize
  // conservativeResize
  // or init at a large size

}

//-------------------------------------------------------------------
// ------------- CONVENIENCE FUNCTIONS -----------------------------
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/*! 
  \brief  Converts from point cloud to NURBS Matrix
  
  \author Benjamin Morrell
  \date 3 April 2018
*/
Matrix_Point3Df Mapping3D::nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){

  // Initialise mesh in format for NURBS
  Matrix_Point3Df mesh(cloud->height,cloud->width);

  for (int i = 0; i < cloud->height; i ++){
    for (int j = 0; j < cloud->width; j++){
      mesh(i,j).x() = cloud->at(j,i).x; // at(col,row)
      mesh(i,j).y() = cloud->at(j,i).y; // at(col,row)
      mesh(i,j).z() = cloud->at(j,i).z; // at(col,row)
    }
  }

  return mesh;
}
Matrix_Point3Df Mapping3D::nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud){

  // Initialise mesh in format for NURBS
  Matrix_Point3Df mesh(cloud->height,cloud->width);

  for (int i = 0; i < cloud->height; i ++){
    for (int j = 0; j < cloud->width; j++){
      mesh(i,j).x() = cloud->at(j,i).x; // at(col,row)
      mesh(i,j).y() = cloud->at(j,i).y; // at(col,row)
      mesh(i,j).z() = cloud->at(j,i).z; // at(col,row)
    }
  }

  return mesh;
}
//-------------------------------------------------------------------
/*! 
  \brief  Converts from point cloud to NURBS Matrix taking only the data indicated by dataIndices
  
  \param cloud  the point cloud
  \param dataIndices indices for the data to take out of the cloud, [starti, startj; endi, endj]

  \author Benjamin Morrell
  \date 3 April 2018
*/
Matrix_Point3Df Mapping3D::nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int, 2, 2>& dataIndices){

  // Size of output matrix
  int ms = dataIndices(1,0) - dataIndices(0,0) + 1;
  int mt = dataIndices(1,1) - dataIndices(0,1) + 1;

  // Initialise mesh in format for NURBS
  Matrix_Point3Df mesh(ms,mt);

  int ii;
  int jj;

  for (int i = 0; i < ms; i ++){
    for (int j = 0; j < mt; j++){
      // Compute indices to get data from cloud
      ii = dataIndices(0,0) + i;
      jj = dataIndices(0,1) + j;
      // Copy data
      mesh(i,j).x() = cloud->at(jj,ii).x; // at(col,row)
      mesh(i,j).y() = cloud->at(jj,ii).y; // at(col,row)
      mesh(i,j).z() = cloud->at(jj,ii).z; // at(col,row)
    }
  }

  return mesh;
}

//-------------------------------------------------------------------
/*! 
  \brief  Converts from point cloud to NURBS Matrix taking only the data indicated by newRowIndices and newColIndices
  
  \param cloud  the point cloud
  \param newRowIndices indices for the rows to take out of the cloud, for each row or col: [start, end]
  \param newColIndices indices for the columns to take out of the cloud, for each row or col: [start, end]

  \author Benjamin Morrell
  \date 3 April 2018
*/
Matrix_Point3Df Mapping3D::nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int, Eigen::Dynamic, 2>& newRowIndices, Eigen::Array<int, Eigen::Dynamic, 2>& newColIndices){

  int nRowOrCol = newRowIndices.rows(); // NUmber of points with indices
  int nPoints;
  int ms;
  int mt;
  int nMeshi;
  bool expandRowsNotCols;
  bool bCountBack;
  int lastVal = -1;

  // Find required size of mesh - smallest size
  if (newRowIndices(0,0) == newRowIndices(0,1)){
    // Constant row - expanding rows
    expandRowsNotCols = true;
    nPoints = (newColIndices.col(1) - newColIndices.col(0)).minCoeff()+1;// Minimum difference between start and end

    // Get number of unique points - the number of rows to take
    // nRowOrCol 

    // Flag for whether to count backward with mesh_i
    if (newRowIndices(0,0)>newRowIndices(nRowOrCol-1,0)){
      bCountBack = true;
    }else{
      bCountBack = false;
    }

    ms = getNumberOfUniquePoints(newRowIndices);
    mt = nPoints;
    nMeshi = ms;

  }else if (newColIndices(0,0) == newColIndices(0,1)) {
    // Constant col - expanding cols
    expandRowsNotCols = false;
    nPoints = (newRowIndices.col(1) - newRowIndices.col(0)).minCoeff()+1;// Minimum difference between start and end
    

    // Get number of unique points - the number of cols to take
    // nRowOrCol = getNumberOfUniquePoints(newColIndices);

    // Flag for whether to count backward with mesh_i
    if (newColIndices(0,0)>newColIndices(nRowOrCol-1,0)){
      bCountBack = true;
    }else{
      bCountBack = false;
    }

    ms = nPoints;
    mt = getNumberOfUniquePoints(newColIndices);
    nMeshi = mt;

  }else{
    cout << "Error in indices, need row or column to be fixed" << endl;
    Matrix_Point3Df mesh(1,1);
    return mesh;
  }

  cout << "Number of rows/cols: " << nRowOrCol << endl;
  cout << "Expanding row is " << expandRowsNotCols << endl;
  cout << "Number of points along row/col is " << nPoints << endl;

  cout << "(m_s, m_t) = (" << ms << ", " << mt << ")\n";

 
  // Initialise mesh in format for NURBS
  Matrix_Point3Df mesh(ms,mt);

  int ii;
  int jj;
  int mesh_i;
  if (bCountBack){
    mesh_i = nMeshi - 1;
  }else{
    mesh_i = 0;
  }

  float stepS;
  float stepT;

  bool newPoint = false;

  // Loop through each row or col
  for (int i = 0; i < nRowOrCol; i++){
    // Step will be zero for S if going along a row, and zero for T if going along a column
    // Step is the number of points from start to end divided by the number of points desired
    stepS = (float)(newRowIndices(i,1) - newRowIndices(i,0))/(float)(nPoints-1); // Should floor the division, so it won't exceed the limit
    stepT = (float)(newColIndices(i,1) - newColIndices(i,0))/(float)(nPoints-1); // Should floor the division, so it won't exceed the limit

    if (expandRowsNotCols){
      // Constant row
      if (newRowIndices(i,0) != lastVal){
        // New Val
        newPoint = true;
        lastVal = newRowIndices(i,0);
      }
    }else{
      // Constant col
      if (newColIndices(i,0) != lastVal){
        // New Val
        newPoint = true;
        lastVal = newColIndices(i,0);
      }
    }

    // cout << "For i = " << i << ", stepS: " << stepS << ", stepT: " << stepT << endl;
    if (newPoint){
      for (int j = 0; j < nPoints; j++ ){
        // Compute indices to get data from cloud
        ii = newRowIndices(i,0) + (int)(j*stepS);
        jj = newColIndices(i,0) + (int)(j*stepT);
        // cout << "(ii,jj) = (" << ii << ", " << jj << ")\n";
        // Copy data
        if (expandRowsNotCols){
          // i is iterating through rows. j along rows
          mesh(mesh_i,j).x() = cloud->at(jj,ii).x; // at(col,row)
          mesh(mesh_i,j).y() = cloud->at(jj,ii).y; // at(col,row)
          mesh(mesh_i,j).z() = cloud->at(jj,ii).z; // at(col,row)
        }else{
          // i is iterating through columns, j along columns
          mesh(j,mesh_i).x() = cloud->at(jj,ii).x; // at(col,row)
          mesh(j,mesh_i).y() = cloud->at(jj,ii).y; // at(col,row)
          mesh(j,mesh_i).z() = cloud->at(jj,ii).z; // at(col,row)
        }
        if (!pcl::isFinite(cloud->at(jj,ii))){
          cout << "\n\nNon-finite mesh point!!\n\n";
        }
      }
      cout << "Mesh_i is: " << mesh_i << endl;
      if (bCountBack){
        mesh_i -= 1;
      }else{
        mesh_i++;
      }
    }
    newPoint = false;
  }

  return mesh;
}

// Function as a helper to get the number of unique points 
int Mapping3D::getNumberOfUniquePoints(Eigen::Array<int, Eigen::Dynamic, 2>& Indices){

  int count = 1;
  int current = Indices(0,0);

  for (int i = 1; i < Indices.rows(); i ++){
    if (current != Indices(i,0)){
      count++;
      current = Indices(i,0);
    }
  }

  cout << "Number of unique points: " << count << endl;

  return count;
}

//-------------------------------------------------------------------
/*! 
  \brief  Converts from NURBS Matrix to pcl Point Cloud
  
  \author Benjamin Morrell
  \date 3 April 2018
*/
void Mapping3D::pointCloudFromNurbsData(Matrix_Point3Df& data, pcl::PointCloud<pcl::PointNormal>::Ptr cloud){

  // pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>(data.cols(),data.rows(),pcl::PointNormal()));

  for (int i = 0; i < data.rows(); i ++){
    for (int j = 0; j < data.cols(); j++){
      cloud->at(j,i).x = data(i,j).x(); // at(col,row)
      cloud->at(j,i).y = data(i,j).y(); // at(col,row)
      cloud->at(j,i).z = data(i,j).z(); // at(col,row)
    }
  }
}

/*! 
  \brief  Converts from NURBS Object to pcl Point Cloud

  \param objID  the id of the object in the map to get the data from
  \param ms     the number of points to create in the s parametric direction
  \param mt     the number of points to create in the t parametric direction
  \param[out]   the output cloud
  
  \author Benjamin Morrell
  \date 4 April 2018
*/
void Mapping3D::pointCloudFromObject3D(int objID, int ms, int mt, pcl::PointCloud<pcl::PointNormal>::Ptr cloud){

  // Get the data points
  objectMap[objID].getSurfacePointCloud(cloud, ms,  mt);

  // Matrix_Point3Df data = objectMap[objID].getSurfacePoints(ms, mt); 

  // cout << "Size of data in pointCloudFromObject3D: (" << data.rows() << ", " << data.cols() << ")\n";

  // cout << "Data (0,0) = " << data(0,0) << endl;

  // // Copy the data
  // for (int i = 0; i < data.rows(); i ++){
  //   for (int j = 0; j < data.cols(); j++){
  //     cloud->at(j,i).x = data(i,j).x(); // at(col,row)
  //     cloud->at(j,i).y = data(i,j).y(); // at(col,row)
  //     cloud->at(j,i).z = data(i,j).z(); // at(col,row)
  //     // try {
  //     //   if (!pcl::isFinite(cloud->at(j,i))){
  //     //     throw 1;
  //     //   }
  //     // } catch(int e){
  //     //   cout << "Non finite point taken from data in pointCloudFromObject3D" << endl;
  //     // }
  //   }
  // }
}


//-------------------------------------------------------------------
/*! 
  \brief  Extracts a row of a matrix

  \param data   from which to extract
  \param row_id row to extract from 
  
  \author Benjamin Morrell
  \date 28 March 2018
*/
Vector_HPoint3Df Mapping3D::getMatRow(Matrix_HPoint3Df data,int row_id){
  // cout << "in getMatRow" << endl;
  int n_ctrl = data.cols();
  
  Vector_HPoint3Df ctrlRow(n_ctrl);
  
  for (int j = 0; j < n_ctrl; j++){
    ctrlRow[j] = data(row_id,j);
  }

  // cout << "end of getMatRow. CtrlRow = " << ctrlRow << endl;
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

//-------------------------------------------------------------------
// ------------- OUTPUT FUNCTIONS -----------------------------
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/*! 
  \brief  Writes object point cloud to file
  
  \author Benjamin Morrell
  \date 3 April 2018
*/
void Mapping3D::writeObjectPCDFile(const char* filename, const int objID, int ms, int mt){

  // Load defaults if no input (ms, mt  = -1 )
  if (ms < 0){
    ms = msSurf;
  }
  if (mt < 0) {
    mt = mtSurf;
  }

  // Get point cloud from NURBS object
  pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC(new pcl::PointCloud<pcl::PointNormal>(mt, ms, pcl::PointNormal()));
  pointCloudFromObject3D(objID, ms, mt, mapObjPC);

  
  pcl::PCDWriter writer;
  writer.write<pcl::PointNormal> (filename, *mapObjPC, false);
}

/*! 
  \brief Fills a ros message with the object of the given object ID
  
  \author Benjamin Morrell
  \date 30 April 2018
*/
void Mapping3D::fillObject3DMessage(int objID, cans_msgs::Object3D& msg){
  
  msg.ID = objID;

  // NURBS parameters
  msg.degU = objectMap[objID].degreeU();
  msg.degV = objectMap[objID].degreeV();
  msg.nCtrlS = objectMap[objID].ctrlPnts().rows();
  msg.nCtrlT = objectMap[objID].ctrlPnts().cols();

  // Knot Vector
  std::vector<float> knotU;
  std::vector<float> knotV;
  for (int i = 0; i < objectMap[objID].knotU().size(); i++){
    knotU.push_back(objectMap[objID].knotU()[i]);
  }
  for (int i = 0; i < objectMap[objID].knotV().size(); i++){
    knotV.push_back(objectMap[objID].knotV()[i]);
  }
  
  msg.knotU = knotU;
  msg.knotV = knotV;

  // Control points
  std::vector<float> controlX;
  std::vector<float> controlY;
  std::vector<float> controlZ;

  for (int i = 0; i < msg.nCtrlS; i++){
    for (int j = 0; j < msg.nCtrlT; j++){
      controlX.push_back(objectMap[objID].ctrlPnts()(i,j).x());
      controlY.push_back(objectMap[objID].ctrlPnts()(i,j).y());
      controlZ.push_back(objectMap[objID].ctrlPnts()(i,j).z());
    }
  }

  msg.ControlPoints_x = controlX;
  msg.ControlPoints_y = controlY;
  msg.ControlPoints_z = controlZ;

}
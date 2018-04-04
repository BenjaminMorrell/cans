
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
    searchThresh(7), numRowsDesired(45), numColsDesired(45),
    maxNanAllowed(10), removeNanBuffer(0), numberOfMetrics(7)
{
  searchThresh[0] = 0.75;
  searchThresh[1] = 0.75;
  searchThresh[2] = 0.75;
  for (int i = 3; i < 7; i++){
    searchThresh[i] = 0.0;  
  }
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

//-------------------------------------------------------------------
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

    Object3D* srfOut;

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

    // Different cases for rows or columns for srf1
    if (extendDirection[0] == 'L' || extendDirection[0] == 'R'){
        // -----------------------------------------------------------------
        // Row-wise expansion
        // -----------------------------------------------------------------
        // Input checks
        if (n_ctrlCheck2 != srf1.ctrlPnts().rows()){
          cout << "Error!: Need to have matched numbers of control points along the join direction\nReturning an empty surface" << endl;
          srfOut = new Object3D();
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
        srfOut = new Object3D(deg,deg,knots1,knotsOut,ctrlOut);
        
    }else{
        // -----------------------------------------------------------------
        // Column-wise expansion
        // -----------------------------------------------------------------
        // Input checks
        if (n_ctrlCheck2 != srf1.ctrlPnts().cols()){
          cout << "Error!: Need to have matched numbers of control points along the join direction\nReturning an empty surface" << endl;
          srfOut = new Object3D();
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
        srfOut = new Object3D(deg,deg,knotsOut,knots1,ctrlOut);

    }
    

    return *srfOut;

}

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
  if (nCtrlNew[0] == 2){
    if ((float)data.rows()/(float)data.cols()*srf.ctrlPnts().rows() > 1.8){
      nCtrlNew[0] = 3;
    }
  }else if (nCtrlNew[1] == 2){
    if ((float)data.cols()/(float)data.rows()*srf.ctrlPnts().cols() > 1.8){
      nCtrlNew[1] = 3;
    }
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
          cloudOut->at(jj,ii).x = cloudIn->at(j,i).x;
          cloudOut->at(jj,ii).y = cloudIn->at(j,i).y;
          cloudOut->at(jj,ii).z = cloudIn->at(j,i).z;

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
  bool noNans = false;
  while (!noNans){
    noNans = averageOutNans(cloudOut, nanIndices);
    cout << "Iteration for averaging Nans" << endl;
  }
  
  
}

/*! 
  \brief  Searches through an organised pointcloud and makes a boolean array with 1s where the point cloud is nan

  \param nanArray   the output boolean array
  \param cloud      a pointer to the point cloud to process
  
  
  \author Benjamin Morrell
  \date 30 March 2018
*/
void Mapping3D::getNanMatrixFromPointCloud(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray, pcl::PointCloud<pcl::PointNormal>::Ptr cloud){

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
bool Mapping3D::averageOutNans(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic>& nanIndices){
  // cout << "Number of cols: " << nanIndices.cols() << endl;

  bool noNans = true;
  cout << "number of nans: " << nanIndices.cols() << endl;

  for (int i = 0; i < nanIndices.cols(); i++){
    if (nanIndices(0,i) == -1){
      // Flag meaning that those terms are not valid (not used)
      break;
    }else if (nanIndices(0,i) == -2){
      // Flag meaning that those terms are not valid (not used), but are not at the end
      continue;
    }
    cout << "Before: " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
    regionAverage(cloud,nanIndices(0,i),nanIndices(1,i));
    cout << "After:  " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
    if (!pcl::isFinite(cloud->at(nanIndices(1,i),nanIndices(0,i)) )){
      // Still nan, keep as flag
      noNans = false;
    }else{
      // Indicate they are not active
      nanIndices(0,i) = nanIndices(1,i) = -2;
    }
  }

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

//-------------------------------------------------------------------
// ------------- UPDATE OBJECT  -----------------------------
//-------------------------------------------------------------------
void Mapping3D::addObject(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, std::vector<float> searchMetrics){

  // Convert data to nurbs data
  Matrix_Point3Df mesh = nurbsDataFromPointCloud(cloud);

  // Add NURBS Object (using default settings)
  objectMap.push_back(Object3D(mesh));

  searchMetrics[0] = objectMap[objectMap.size()-1].getCentre().x();
  searchMetrics[1] = objectMap[objectMap.size()-1].getCentre().y();
  searchMetrics[2] = objectMap[objectMap.size()-1].getCentre().z();
  searchMetrics[3] = objectMap[objectMap.size()-1].getObjSize();

  // check search metrics size

  // Add search metrics
  objectMetrics.push_back(searchMetrics);

}

void Mapping3D::updateObject(Object3D& mapObj, pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC){
  // TDBM - function to get point cloud from NURBS
  // TDBM just input an ID and get the data from the object map
  
  // Split new surface observation
  SplitSurface ss;

  ss.splitNewSurfaceObservation(mapObjPC, obsObjPC);

  cout << "Extend Direction is: " << ss.extendDirection << endl;

  // Get matrix from new data
  Matrix_Point3Df mesh = nurbsDataFromPointCloud(obsObjPC, ss.newDataIndices);

  cout << "Mesh2 size is (" << mesh.rows() << ", " << mesh.cols() << ")\n";

  std::vector<int> nCtrlNew = computeNumberOfControlPoints(ss.extendDirection, mesh, mapObj);

  cout << "new control points: " << nCtrlNew[0] << ", " << nCtrlNew[1] << endl;

  // Create Object for new surface? (nknots = deg + nctrl + 1)
  Object3D obj(mesh, order[0], order[1], nCtrlNew[0], nCtrlNew[1]);

  // Join Surfaces
  Object3D objJoin = joinSurfaces(mapObj,obj,ss.extendDirection);

  cout << "Point (0.3,0.3) on obj3: " << objJoin.pointAt(0.3,0.3) << endl;

  objJoin.writeVRML("mesh3Test.wrl",Color(255,100,255),50,80);

}

//-------------------------------------------------------------------
// ------------- DATA ASSOCIATION  -----------------------------
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/*! 
  \brief  Converts from point cloud to NURBS Matrix

  uses: objectMetrics
        numberOfMetrics
        searchThresh
  
  \author Benjamin Morrell
  \date 3 April 2018
*/
int Mapping3D::dataAssociation(std::vector<float> searchMetrics){

  // Input
  // Eigen::Array<float, 1, Eigen::Dynamic> searchMetrics(1,7);

  // From object map [x, y, z, size, r, g, b]. Init here with 10 objects as an example
  // Eigen::Array<float,Eigen::Dynamic, Eigen::Dynamic> objectMetrics(10,7);

  // std::vector<float> searchThresh;

  int objID;
  float dist;

  // Distances array
  Eigen::Array<float, 1, Eigen::Dynamic> distances(1,numberOfMetrics);
  distances.setZero(numberOfMetrics);

  Eigen::Array<bool, 1, Eigen::Dynamic> activeObjects(1,numberOfMetrics);
  activeObjects.setOnes(numberOfMetrics);

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
    while (activeObjects.count() > 1){
      // Get index of the first remaining match
      activeObjects.maxCoeff(&i);
      
      // Get the distance 
      dist = distances[i];

      if (dist < minDist){
        // Select objID to get min
        objID = i;
        // update min dist
        minDist = dist;
      }
      // Update activeObjects to get next search
      activeObjects[i] = false;
    }
  }

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
Matrix_Point3Df Mapping3D::nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int, 2, 2> dataIndices){

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
  \brief  Converts from NURBS Matrix to pcl Point Cloud
  
  \author Benjamin Morrell
  \date 3 April 2018
*/
pcl::PointCloud<pcl::PointNormal> Mapping3D::pointCloudFromNurbsData(Matrix_Point3Df data){

  pcl::PointCloud<pcl::PointNormal> cloud(data.cols(),data.rows(),pcl::PointNormal());

  for (int i = 0; i < data.rows(); i ++){
    for (int j = 0; j < data.cols(); j++){
      cloud.at(j,i).x = data(i,j).x(); // at(col,row)
      cloud.at(j,i).y = data(i,j).y(); // at(col,row)
      cloud.at(j,i).z = data(i,j).z(); // at(col,row)
    }
  }

  return cloud;// TDBM think of how to do this more efficiently with pointers - depends on how it is used...
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
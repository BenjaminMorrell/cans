
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
    search_thresh {0.75,0.75,0.75,0,0,0,0,0,0}
{
}


//-------------------------------------------------------------------
/*! 
  \brief  Find the centre of data, to be used for data association

  \author Benjamin Morrell
  \date 23 March 2018
*/
Point3Df Mapping3D::compute_centre_of_data(const Matrix_Point3Df &scan){

    Point3Df centre(0.0,0.0,0.0);

    // Define step size to only use a subset of points for sampling the mean
    int step_size(1);

    // Add up all the points
    for (int i = 0; i < scan.rows(); i++){
        for (int j = 0; j < scan.cols(); j++){
            // Compute the mean for x, y and z
            centre += scan(i,j);
        }
    }
    // Divide by the number of terms to get the average
    centre /= (float)(scan.rows()*scan.cols());
    // centre /= static_cast<float>(scan.rows()*scan.cols())



    // return by value as it is a local variable, and not very large
    return centre;
}

//-------------------------------------------------------------------
/*! 
  \brief  Form an ordered mesh from input scan data

  \author Benjamin Morrell
  \date 23 March 2018
*/
Matrix_Point3Df Mapping3D::mesh_from_scan(const Matrix_Point3Df scan){

    

    return scan;

    // Try to pass back a pointer to the data?

}




//-------------------------------------------------------------------
// PRIVATE METHODS
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
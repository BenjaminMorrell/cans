#include "cans/object3D.h"

using namespace PLib ; 
//-------------------------------------------------------------------
/*! 
  \brief  Default Constructor

  \author Benjamin Morrell
  \date 23 March 2018
*/

Object3D::Object3D():
  centre(0,0,0), 
  color(0,0,0),
  obj_size(0.0),
  GIS_texture(0.0),
  GIS_roughness(0.0)
{}


// //-------------------------------------------------------------------
// /*! 
//   \brief  General constructor with default inputs

//   \warning  Could this be used as the only constructor? 

//   \author Benjamin Morrell
//   \date 23 March 2018
// */
// Object3D::Object3D(PlNurbsSurfacef &nS, Point3Df &cp, Point3Df &col, float obj_size, float GIS1, float GIS2)  :
//   PlNurbsSurfacef(nS),
//   centre(cp),
//   obj_size(obj_size),
//   color(col),
//   GIS_texture(GIS1),
//   GIS_roughness(GIS2) {}



//-------------------------------------------------------------------
/*! 
  \brief  Least Squares fit constructor

  \param   Q  a matrix of 3D points
  \param  pU  the degree of interpolation in the U direction
  \param  pV  the degree of interpolation in the V direction
  \param  nU  the number of points in the U direction
  \param  nV  the number of poitns in the V direction

  \author Benjamin Morrell
  \date 23 March 2018
*/
Object3D::Object3D(const Matrix_Point3Df& Q, int pU, int pV, int nU, int nV){
  // do a least sqaures fit to the input data
  leastSquares(Q, pU, pV, nU, nV);

  // Compute the centre from the data (or the control points?)
  computeCentreFromData(Q); // sets the centre variable
  // computeCentreFromControlPoints(); // sets the centre variable

  // Compute the size from the data
  computeSizeFromData(Q); // sets the object_size variable
  // computeSizeFromControlPoints(); // sets the object_size variable

}

//-------------------------------------------------------------------
/*! 
  \brief  Least Squares fit constructor

  \param   Q  a matrix of 3D points

  \warning Uses default values for nurbs parameters

  \author Benjamin Morrell
  \date 23 March 2018
*/
Object3D::Object3D(const Matrix_Point3Df& Q){
  // do a least sqaures fit to the input data
  // Default order is 3
  // Default number of knots is 10
  leastSquares(Q, 3, 3, 10, 10); // TODO: Make the defaults adapt with the size of the matrix

  // Compute the centre from the data (or the control points?)
  computeCentreFromData(Q); // sets the centre variable
  // computeCentreFromControlPoints(); // sets the centre variable

  // Compute the size from the data
  computeSizeFromData(Q); // sets the object_size variable
  // computeSizeFromControlPoints(); // sets the object_size variable

}



//-------------------------------------------------------------------
/*! 
  \brief  Computing the centre of a matrix of input data points

  \param   Q  a matrix of 3D points

  \author Benjamin Morrell
  \date 23 March 2018
*/
void Object3D::computeCentreFromData(const Matrix_Point3Df& scan, int step_size){
  // Need to check that you can edit a point like this
  centre.data[0] = 0.0;
  centre.data[1] = 0.0;
  centre.data[2] = 0.0;

  // Add up all the points
  for (int i = 0; i < scan.rows(); i = i + step_size){
    for (int j = 0; j < scan.cols(); j = j + step_size){
      // Compute the mean for x, y and z
      centre += scan(i,j);
    }
  }
  // Divide by the number of terms to get the average
  centre /= (float)(scan.rows()*scan.cols());
}

//-------------------------------------------------------------------
/*! 
  \brief  Computing the centre from control points

  \param   Q  a matrix of 3D points

  \author Benjamin Morrell
  \date 23 March 2018
*/
void Object3D::computeCentreFromControlPoints(){
  // Need to check that you can edit a point like this
  centre.data[0] = 0.0;
  centre.data[1] = 0.0;
  centre.data[2] = 0.0;

  Matrix_HPoint3Df control_points = this->ctrlPnts();
  
  // Add up all the points
  for (int i = 0; i < control_points.rows(); i++){
    for (int j = 0; j < control_points.cols(); j++){
      // Compute the mean for x, y and z
      centre.x() += control_points(i,j).x();
      centre.y() += control_points(i,j).y();
      centre.z() += control_points(i,j).z();
    }
  }
  // Divide by the number of terms to get the average
  centre /= (float)(control_points.rows()*control_points.cols());
}

//-------------------------------------------------------------------
/*! 
  \brief  Computing the size from input data

  \author Benjamin Morrell
  \date 23 March 2018
*/
void Object3D::computeSizeFromData(const Matrix_Point3Df& Q){

  Point3Df minVals(999999,999999,999999);
  Point3Df maxVals(-999999,-999999,-999999);

  // Add up all the points
  for (int i = 0; i < Q.rows(); i++){
    for (int j = 0; j < Q.cols(); j++){
      if (Q(i,j).x() < minVals.x()){
        minVals.x() = Q(i,j).x();
      }
      if (Q(i,j).x() > maxVals.x()){
        maxVals.x() = Q(i,j).x();
      }
      if (Q(i,j).y() < minVals.y()){
        minVals.y() = Q(i,j).y();
      }
      if (Q(i,j).y() > maxVals.y()){
        maxVals.y() = Q(i,j).y();
      }
      if (Q(i,j).z() < minVals.z()){
        minVals.z() = Q(i,j).z();
      }
      if (Q(i,j).z() > maxVals.z()){
        maxVals.z() = Q(i,j).z();
      }
    }
  }

  // Max delta
  Point3Df delta = maxVals - minVals;

  // Size as the maximum dimenation TODO: Come up with a better metric
  obj_size = std::max(delta.x(),std::max(delta.y(),delta.z()));
  

}

//-------------------------------------------------------------------
/*! 
  \brief  Computing the size from control points

  \author Benjamin Morrell
  \date 23 March 2018
*/
void Object3D::computeSizeFromControlPoints(){

  Point3Df minVals(999999,999999,999999);
  Point3Df maxVals(-999999,-999999,-999999);

  Matrix_HPoint3Df control_points = this->ctrlPnts();

  // Add up all the points
  for (int i = 0; i < control_points.rows(); i++){
    for (int j = 0; j < control_points.cols(); j++){
      if (control_points(i,j).x() < minVals.x()){
        minVals.x() = control_points(i,j).x();
      }
      if (control_points(i,j).x() > maxVals.x()){
        maxVals.x() = control_points(i,j).x();
      }
      if (control_points(i,j).y() < minVals.y()){
        minVals.y() = control_points(i,j).y();
      }
      if (control_points(i,j).y() > maxVals.y()){
        maxVals.y() = control_points(i,j).y();
      }
      if (control_points(i,j).z() < minVals.z()){
        minVals.z() = control_points(i,j).z();
      }
      if (control_points(i,j).z() > maxVals.z()){
        maxVals.z() = control_points(i,j).z();
      }
    }
  }

  // Max delta
  Point3Df delta = maxVals - minVals;

  // Size as the maximum dimenation TODO: Come up with a better metric
  obj_size = std::max(delta.x(),std::max(delta.y(),delta.z()));

}


//-------------------------------------------------------------------
/*! 
  \brief  Get the centre

  \author Benjamin Morrell
  \date 23 March 2018
*/
Point3Df& Object3D::getCentre(){

  // do I need to do this->centre?
  return centre;
}

//-------------------------------------------------------------------
/*! 
  \brief  Get the color

  \author Benjamin Morrell
  \date 23 March 2018
*/
Point3Df& Object3D::getColor(){

  // do I need to do this->centre?
  return color;
}

//-------------------------------------------------------------------
/*! 
  \brief  Get the object size parameter

  \author Benjamin Morrell
  \date 23 March 2018
*/
float Object3D::getObjSize(){

  // do I need to do this->centre?
  return obj_size;
}





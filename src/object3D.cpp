#include "cans/object3D.h"

using namespace PLib ; 
//-------------------------------------------------------------------
/*! 
  \brief  Default Constructor

  \author Benjamin Morrell
  \date 23 March 2018
*/

Object3D::Object3D():
  NurbsSurfacef(),
  centre(0,0,0), 
  color(0,0,0),
  objSize(0.0),
  GISTexture(0.0),
  GISRoughness(0.0)
{}

//-------------------------------------------------------------------
/*! 
  \brief  Copy Constructor

  \author Benjamin Morrell
  \date 23 March 2018
*/

Object3D::Object3D(const Object3D & obj):
  NurbsSurfacef(obj),
  centre(obj.centre), 
  color(obj.color),
  objSize(obj.objSize),
  GISTexture(obj.GISTexture),
  GISRoughness(obj.GISRoughness)
{}


// //-------------------------------------------------------------------
// /*! 
//   \brief  General constructor with default inputs

//   \warning  Could this be used as the only constructor? 

//   \author Benjamin Morrell
//   \date 23 March 2018
// */
// Object3D::Object3D(PlNurbsSurfacef &nS, Point3Df &cp, Point3Df &col, float objSize, float GIS1, float GIS2)  :
//   PlNurbsSurfacef(nS),
//   centre(cp),
//   objSize(objSize),
//   color(col),
//   GISTexture(GIS1),
//   GISRoughness(GIS2) {}



//-------------------------------------------------------------------
/*! 
  \brief  Least Squares fit constructor

  \param   Q  a matrix of 3D points
  \param  pU  the degree of interpolation in the U direction
  \param  pV  the degree of interpolation in the V direction
  \param  nU  the number of control points in the U direction
  \param  nV  the number of control points in the V direction

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
  \brief  NURBS initialised constructor

  \param   Q  a matrix of 3D points

  \warning Uses default values for nurbs parameters

  \author Benjamin Morrell
  \date 23 March 2018
*/
Object3D::Object3D(int pU, int pV, Vector_FLOAT& UVec, Vector_FLOAT& Vvec, Matrix_HPoint3Df& ctrlPnts) :
  NurbsSurfacef(pU, pV, UVec, Vvec, ctrlPnts)
{

  // Compute the centre from the data (or the control points?)
  // computeCentreFromData(Q); // sets the centre variable
  computeCentreFromControlPoints(); // sets the centre variable

  // Compute the size from the data
  // computeSizeFromData(Q); // sets the object_size variable
  computeSizeFromControlPoints(); // sets the object_size variable

}

Object3D::~Object3D(){
  ;
}

//-------------------------------------------------------------------
/*! 
  \brief  Assignment operator 
  
  TDBM - this may not copy the nurbsSruface properties...

  \author Benjamin Morrell
  \date 23 March 2018
*/
void Object3D::operator = (const Object3D& obj){
  NurbsSurfacef::operator = (obj);
  centre = obj.centre; 
  color = obj.color;
  objSize = obj.objSize;
  GISTexture = obj.GISTexture;
  GISRoughness = obj.GISRoughness;
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
  objSize = std::max(delta.x(),std::max(delta.y(),delta.z()));
  

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
  objSize = std::max(delta.x(),std::max(delta.y(),delta.z()));

}

//-------------------------------------------------------------------
/*! 
  \brief  Load a nurbs and compute the centre etc. 

  \warning color will be zeros - have not copied them

  \author Benjamin Morrell
  \date 17 April 2018
*/
void Object3D::readObject3D(const char* filename){
  // Read the file
  read(filename);

  // Compute the centre 
  computeCentreFromControlPoints();

  // Compute the size
  computeSizeFromControlPoints();

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
  return objSize;
}

//-------------------------------------------------------------------
/*! 
  \brief  Get data from point cloud

  \param ms number of points in s parametric direction
  \param mt number of points in t parametric direction

  \author Benjamin Morrell
  \date 23 March 2018
*/
Matrix_Point3Df Object3D::getSurfacePoints(int ms, int mt){

  Matrix_Point3Df data(ms,mt);

  // cout << "ms is " << ms << ", mt is " << mt << ", size of data is (" << data.rows() << ", " << data.cols() << ")\n";

  // Create parameter vectors
  Vector_FLOAT paramsS(ms);
  float val(0.0001);
  float step = (1.0-val)/ms;

  for (int k = 0; k < ms; k++){
    paramsS[k] = val;
    val += step;
  }

  Vector_FLOAT paramsT(mt);
  val = 0.0001;
  step = (1.0-val)/mt;

  for (int k = 0; k < mt; k++){
    paramsT[k] = val;
    val += step;
  }

  // cout << "params are: " << paramsS << "\n\n" << paramsT << endl;

  // Fill matrix
  for (int i = 0; i < ms; i++){
    for (int j = 0; j < mt; j++){
      // Point cloud with normals... 
      data(i,j).x() = this->pointAt(paramsS[i],paramsT[j]).x();
      data(i,j).y() = this->pointAt(paramsS[i],paramsT[j]).y();
      data(i,j).z() = this->pointAt(paramsS[i],paramsT[j]).z();
      // cout << "Data from NURBS object is: " << data(i,j).x() << endl;
    }
  }
  return data;
}

// //-------------------------------------------------------------------
// /*! 
//   \brief  Get data from point cloud

//   \author Benjamin Morrell
//   \date 23 March 2018
// */
// pcl::PointCloud<pcl::PointNormal>::Ptr Object3D::getSurfacePointCloud(int ms = 45, int mt = 45){
//   pcl::PointCloud<pcl::PointNormal>::Ptr new_cloud( new pcl::PointCloud<pcl::PointNormal>(ms,mt,PointNormal()));
//   for (int i = 0; i< ms; i++){
//     for (int j = 0; j < mt; j++){
//       // Point cloud with normals... 
//       new_cloud_n->at(j,i).x = this->pointAt(params[i],params[j]).x();
//       new_cloud_n->at(j,i).y = this->pointAt(params[i],params[j]).y();
//       new_cloud_n->at(j,i).z = this->pointAt(params[i],params[j]).z();

//       new_cloud_n->at(j,i).normal_x = this->normal(params[i],params[j]).x();
//       new_cloud_n->at(j,i).normal_y = this->normal(params[i],params[j]).y();
//       new_cloud_n->at(j,i).normal_z = this->normal(params[i],params[j]).z();
//     }
//   }
// }





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
Object3D::Object3D(int pU, int pV, Vector_FLOAT& Uvec, Vector_FLOAT& Vvec, Matrix_HPoint3Df& ctrlPnts) :
  NurbsSurfacef(pU, pV, Uvec, Vvec, ctrlPnts)
{

  // Compute the centre from the data (or the control points?)
  // computeCentreFromData(Q); // sets the centre variable
  computeCentreFromControlPoints(); // sets the centre variable

  // Compute the size from the data
  // computeSizeFromData(Q); // sets the object_size variable
  computeSizeFromControlPoints(); // sets the object_size variable

}


Object3D::Object3D(int pU, int pV, Eigen::Array<float,1,Eigen::Dynamic>& Uvec, Eigen::Array<float,1,Eigen::Dynamic>& Vvec, Eigen::Array<float,3,Eigen::Dynamic>& ctrlPnts, int nCtrlS, int nCtrlT):
  NurbsSurfacef(),
  centre(0,0,0), 
  color(0,0,0),
  objSize(0.0),
  GISTexture(0.0),
  GISRoughness(0.0)
{

  if (Uvec.size() != nCtrlS + pU + 1){
    cout << "Error, knot vector and control point vector mis-match" << endl;
    return;
  }else if (Vvec.size() != nCtrlT + pV + 1) {
    cout << "Error, knot vector and control point vector mis-match" << endl;
    return;
  }

  this->resize(nCtrlS,nCtrlT, pU, pV);
  
  // Set up knots
  int nKnotU = Uvec.size();// or may have to do cols()
  int nKnotV = Vvec.size();// or may have to do cols()

  Vector_FLOAT UvecNURBS(nKnotU);
  Vector_FLOAT VvecNURBS(nKnotV);

  for (int i = 0; i < nKnotU; i++){
    UvecNURBS[i] = Uvec(i);
  }
  for (int i = 0; i < nKnotV; i++){
    VvecNURBS[i] = Vvec(i);
  }

  // Set up Control points 
  Matrix_HPoint3Df controlP(nCtrlS,nCtrlT);

  int index = 0;
  for (int i = 0; i < nCtrlS; i++){
    for (int j = 0; j < nCtrlT; j ++){
      controlP(i,j).x() = ctrlPnts(0,index);
      controlP(i,j).y() = ctrlPnts(1,index);
      controlP(i,j).z() = ctrlPnts(2,index);
      controlP(i,j).w() = 1.0f;
      index ++;
    }
  }

  // Set values in the member variables:
  degU = pU;
  degV = pV;
  U = UvecNURBS;
  V = VvecNURBS;
  P = controlP;
  
  // Compute the centre from the control points
  computeCentreFromControlPoints(); // sets the centre variable

  // Compute the size from the control points
  computeSizeFromControlPoints(); // sets the object_size variable

}

//-------------------------------------------------------------------
/*! 
  \brief  destructor

  \author Benjamin Morrell
  \date 23 March 2018
*/
Object3D::~Object3D(){
  ;
}

//-------------------------------------------------------------------
/*! 
  \brief  Update the object - for use as a python bound object (why Eigen is used)
  
  TDBM - this may not copy the nurbsSruface properties...

  \author Benjamin Morrell
  \date 23 March 2018
*/
void Object3D::updateObject3D(int pU, int pV, Eigen::Array<float,1,Eigen::Dynamic>& Uvec, Eigen::Array<float,1,Eigen::Dynamic>& Vvec, Eigen::Array<float,3,Eigen::Dynamic>& ctrlPnts, int nCtrlS, int nCtrlT)
{

  if (Uvec.size() != nCtrlS + pU + 1){
    cout << "Error, knot vector and control point vector mis-match" << endl;
    return;
  }else if (Vvec.size() != nCtrlT + pV + 1) {
    cout << "Error, knot vector and control point vector mis-match" << endl;
    return;
  }

  this->resize(nCtrlS,nCtrlT, pU, pV);

  cout << "resized in update Object3D" << endl;
  
  // Set up knots
  int nKnotU = Uvec.size();// or may have to do cols()
  int nKnotV = Vvec.size();// or may have to do cols()

  Vector_FLOAT UvecNURBS(nKnotU);
  Vector_FLOAT VvecNURBS(nKnotV);

  for (int i = 0; i < nKnotU; i++){
    UvecNURBS[i] = Uvec(i);
  }
  for (int i = 0; i < nKnotV; i++){
    VvecNURBS[i] = Vvec(i);
  }

  // Set up Control points 
  Matrix_HPoint3Df controlP(nCtrlS,nCtrlT);

  int index = 0;
  for (int i = 0; i < nCtrlS; i++){
    for (int j = 0; j < nCtrlT; j ++){
      controlP(i,j).x() = ctrlPnts(0,index);
      controlP(i,j).y() = ctrlPnts(1,index);
      controlP(i,j).z() = ctrlPnts(2,index);
      controlP(i,j).w() = 1.0f;
      index ++;
    }
  }

  // Set values in the member variables:
  degU = pU;
  degV = pV;
  U = UvecNURBS;
  V = VvecNURBS;
  P = controlP;
  
  // Compute the centre from the control points
  computeCentreFromControlPoints(); // sets the centre variable

  // Compute the size from the control points
  computeSizeFromControlPoints(); // sets the object_size variable

}

//-------------------------------------------------------------------
/*! 
  \brief  Update the object - for use as the c++ listener
  
  TDBM - this may not copy the nurbsSruface properties...

  \author Benjamin Morrell
  \date 23 March 2018
*/
void Object3D::updateObject3DCPP(int pU, int pV, std::vector<float>& Uvec, std::vector<float>& Vvec, std::vector<float>& ctrlPntsX, std::vector<float>& ctrlPntsY, std::vector<float>& ctrlPntsZ, int nCtrlS, int nCtrlT)
{

  // if (Uvec.size() != nCtrlS + pU + 1){
  //   cout << "Error, knot vector and control point vector mis-match" << endl;
  //   return;
  // }else if (Vvec.size() != nCtrlT + pV + 1) {
  //   cout << "Error, knot vector and control point vector mis-match" << endl;
  //   return;
  // }

  this->resize(nCtrlS,nCtrlT, pU, pV);
  
  // Set up knots
  int nKnotU = Uvec.size();// or may have to do cols()
  int nKnotV = Vvec.size();// or may have to do cols()

  Vector_FLOAT UvecNURBS(nKnotU);
  Vector_FLOAT VvecNURBS(nKnotV);

  for (int i = 0; i < nKnotU; i++){
    UvecNURBS[i] = Uvec[i];
  }
  for (int i = 0; i < nKnotV; i++){
    VvecNURBS[i] = Vvec[i];
  }

  // Set up Control points 
  Matrix_HPoint3Df controlP(nCtrlS,nCtrlT);

  int index = 0;
  for (int i = 0; i < nCtrlS; i++){
    for (int j = 0; j < nCtrlT; j ++){
      controlP(i,j).x() = ctrlPntsX[index];
      controlP(i,j).y() = ctrlPntsY[index];
      controlP(i,j).z() = ctrlPntsZ[index];
      controlP(i,j).w() = 1.0f;
      index ++;
    }
  }

  // Set values in the member variables:
  degU = pU;
  degV = pV;
  U = UvecNURBS;
  V = VvecNURBS;
  P = controlP;
  
  // Compute the centre from the control points
  computeCentreFromControlPoints(); // sets the centre variable

  // Compute the size from the control points
  computeSizeFromControlPoints(); // sets the object_size variable

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
  \brief  Get data from NURBS surface

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
      // try {
      //   if (!std::isfinite(data(i,j).x())){
      //     throw 1;
      //   }
      // }catch (int e){
      //   cout << "Non finite value in getSurfacePoints, at (i,j) = (" << i << ", " << j << ")." << endl;
      // }
    }
  }
  return data;
}

//-------------------------------------------------------------------
/*! 
  \brief  Get point cloud from NURBS surface

  \author Benjamin Morrell
  \date 19 April 2018
*/
void Object3D::getSurfacePointCloud( pcl::PointCloud<pcl::PointNormal>::Ptr cloud,int ms, int mt){

  // TODO - input check if NURBS Is empty TDBM
  // Assume cloud is initialised at the correct size
  // // Initialise cloud?
  // pcl::PointCloud<pcl::PointNormal> new_cloud( new pcl::PointCloud<pcl::PointNormal>(ms,mt,PointNormal()));
  try{
    if (cloud->height != ms || cloud->width != mt){  
      throw std::string("Cloud Dimensions Wrong");
    }
  }catch (std::string e) {
    cout << "Error, need cloud to be correct dimensions in getSurfacePointCloud" << endl;
    cout << "Cloud dimensions are, h: " << cloud->height << " w: " << cloud->width << endl;
  }

  // Create parameter vectors
  Vector_FLOAT paramsS(ms);
  float val(0.0001);
  float step = (1.0-val*2.0)/(float)(ms-1.0);

  for (int k = 0; k < ms; k++){
    paramsS[k] = val;
    val += step;
  }

  Vector_FLOAT paramsT(mt);
  val = 0.0001;
  step = (1.0-val*2.0)/(float)(mt-1.0);

  for (int k = 0; k < mt; k++){
    paramsT[k] = val;
    val += step;
  }

  // cout << "ParamS: " << paramsS << endl;
  // cout << "ParamT: " << paramsT << endl;

  // Fill the point cloud
  
    for (int i = 0; i< ms; i++){
      for (int j = 0; j < mt; j++){
        // Point cloud with normals... 
        try{
        cloud->at(j,i).x = this->pointAt(paramsS[i],paramsT[j]).x();
        cloud->at(j,i).y = this->pointAt(paramsS[i],paramsT[j]).y();
        cloud->at(j,i).z = this->pointAt(paramsS[i],paramsT[j]).z();
        
        cloud->at(j,i).normal_x = this->normal(paramsS[i],paramsT[j]).x();
        cloud->at(j,i).normal_y = this->normal(paramsS[i],paramsT[j]).y();
        cloud->at(j,i).normal_z = this->normal(paramsS[i],paramsT[j]).z();
        }catch(...){
          try{
            cout << "cloud point is: " << cloud->at(j,i).x << endl;
          }catch (...){
            cout << "cloud point failed" << endl;
          }
          try{
            cout << "surface point is: " << this->pointAt(paramsS[i],paramsT[j]).x() << endl;
          }catch (...){
            cout << "surface point failed" << endl;
            cout << "Params are: " << paramsS[i] << ", " << paramsT[j] << endl;
          }
          // cout << "Filling point failed" << endl;
        }
      }
    }
  // }catch (...) {
  //   cout << "Filling the point cloud failed!" << endl;
  // }
}


//-------------------------------------------------------------------
/*! 
  \brief  Get signed distance from the point to the surface

  \param[in]  query      the points to get the distances for
  \param[out]            the distance to the point
  \param      ms         The number of points to produce from the surface in the s direction
  \param      mt         The number of points to produce from the surface in the t direction


  \author Benjamin Morrell
  \date 19 April 2018
*/
float Object3D::getDistanceFromPointToSurface(Eigen::Vector3f& query, int ms, int mt){
  // TODO - input check if NURBS Is empty TDBM

  // // Settings for points to get. TDBM set this smartly
  // int ms = 25;
  // int mt = 25;

  // Get point cloud from NURBS
  // Initialise point cloud
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud( new pcl::PointCloud<pcl::PointNormal>(mt,ms,pcl::PointNormal()));

  // Fill point cloud from NURBS
  // try{
    getSurfacePointCloud(cloud, ms, mt);
    cout << "ms, mt = " << ms << ", " << mt << endl;
    cout << "point cloud size is " << cloud->size() << endl;
  // }
  // catch (...){
  //   cout << "Get Surface point Cloud Failed" << endl;
  //   cout << "ms, mt = " << ms << ", " << mt << endl;
  //   cout << "point cloud size is " << cloud->size() << ": (" << cloud->height << ", " << cloud->width << ")" << endl;
  //   return -999999.9;
  // }

  // Create point cloud from input point
  pcl::PointCloud<pcl::PointNormal>::Ptr query_cloud( new pcl::PointCloud<pcl::PointNormal>(1,1,pcl::PointNormal()));

  query_cloud->at(0).x = query(0);
  query_cloud->at(0).y = query(1);
  query_cloud->at(0).z = query(2);

  cout << "Query cloud is " << query_cloud->at(0) << endl;
  cout << "size is: " << query_cloud->size() << endl;

  // Set up correspondences
  pcl::registration::CorrespondenceEstimation<pcl::PointNormal, pcl::PointNormal> corrEst;
  corrEst.setInputSource (query_cloud);// Gets a corresponding point for every point in the Source
  corrEst.setInputTarget (cloud); // Target for the source to find corresponding points in
  
  pcl::CorrespondencesPtr all_correspondencesPtr (new pcl::Correspondences);
  pcl::Correspondences& all_correspondences = *all_correspondencesPtr;
  
  cout << "Initialised correspondences" << endl;

  try {
    // Determine all correspondences
    corrEst.determineCorrespondences (*all_correspondencesPtr);
  }
  catch (...) {
    cout << "Correspondence estimation failed " << endl;
    return -999999.9;
  }

  cout << "Correspondences determined " << endl;

  cout << "Correspondence index query: " << all_correspondences[0].index_query << endl;
  cout << "Correspondence match index: " << all_correspondences[0].index_match << endl;
  cout << "Correspondence distance: " << all_correspondences[0].distance << endl;
  cout << "Correspondence weight: " << all_correspondences[0].weight << endl;

  // Determine the sign
  // Get the point of the match and the normal at that point
  Eigen::Vector3f matchPoint;
  Eigen::Vector3f matchNormal;
  matchPoint(0) = cloud->at(all_correspondences[0].index_match).x;
  matchPoint(1) = cloud->at(all_correspondences[0].index_match).y;
  matchPoint(2) = cloud->at(all_correspondences[0].index_match).z;
  matchNormal(0) = -cloud->at(all_correspondences[0].index_match).normal_x;
  matchNormal(1) = -cloud->at(all_correspondences[0].index_match).normal_y;
  matchNormal(2) = -cloud->at(all_correspondences[0].index_match).normal_z;

  cout << "Normal is: " << matchNormal << endl;

  // Get the vector from the surface to the test point 
  Eigen::Vector3f delta = query - matchPoint;

  // Get the dot product - negative if the point is inside the surface (a dot b) = |a||b|cos(theta). |theta| > 90 deg means it is inside surface - results in a negative dot product
  float dotProd = query.dot(matchNormal);
  cout << "Dot product result is: " << dotProd << endl;

  // TDBM buffer on signed check?
  if (dotProd < 0){
    // Point inside object
    cout << "Signed distance is: " << all_correspondences[0].distance*-1.0 << endl;
    return std::sqrt(all_correspondences[0].distance)*-1.0;
  }else{
    // Point Outside object
    cout << "Signed distance is: " << all_correspondences[0].distance << endl;
    return std::sqrt(all_correspondences[0].distance);
  }

}


//-------------------------------------------------------------------
/*! 
  \brief  Get signed distance from a set of the points to the surface

  \param[in]  query      the points to get the distances for
  \param[out] distAndGrad  the distances output and the gradient [dist;gradx;grady;gradz]
  \param      ms         The number of points to produce from the surface in the s direction
  \param      mt         The number of points to produce from the surface in the t direction


  Gradient is approximated by using the vector from the matched surface point to the query point

  \author Benjamin Morrell
  \date 19 April 2018
*/
Eigen::Array<float,4,Eigen::Dynamic> Object3D::getBatchDistanceFromPointsToSurface(Eigen::Array<float,3,Eigen::Dynamic>& query, int ms, int mt){
  // TODO - input check if NURBS Is empty TDBM

  Eigen::Array<float,4,Eigen::Dynamic> distAndGrad(4,query.cols());

  // if (query.cols() != distances.cols()){
  //   cout << "Error: distances input needs to be same size as query" << endl;
  //   return;
  // }

  // // Settings for points to get. TDBM set this smartly
  // int ms = 25;
  // int mt = 25;

  // Get point cloud from NURBS
  // Initialise point cloud
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud( new pcl::PointCloud<pcl::PointNormal>(mt,ms,pcl::PointNormal()));

  // Fill point cloud from NURBS
  getSurfacePointCloud(cloud, ms, mt);

  // Create point cloud from input point
  int nDataPoints = query.cols();
  pcl::PointCloud<pcl::PointNormal>::Ptr query_cloud( new pcl::PointCloud<pcl::PointNormal>(1,nDataPoints,pcl::PointNormal()));

  cout << "Query has " << nDataPoints << " data points." << endl;

  for (int j = 0; j < nDataPoints; j++){
    query_cloud->at(j).x = query(0,j);
    query_cloud->at(j).y = query(1,j);
    query_cloud->at(j).z = query(2,j);
  }
  
  // Set up correspondences
  pcl::registration::CorrespondenceEstimation<pcl::PointNormal, pcl::PointNormal> corrEst;
  corrEst.setInputSource (query_cloud);// Gets a corresponding point for every point in the Source
  corrEst.setInputTarget (cloud); // Target for the source to find corresponding points in
  
  pcl::CorrespondencesPtr all_correspondencesPtr (new pcl::Correspondences);
  pcl::Correspondences& all_correspondences = *all_correspondencesPtr;
  
  cout << "Initialised correspondences" << endl;

  // Determine all correspondences
  corrEst.determineCorrespondences (*all_correspondencesPtr);

  cout << "Correspondences determined " << endl;

  cout << "Correspondence 0 index query: " << all_correspondences[0].index_query << endl;
  cout << "Correspondence 0 match index: " << all_correspondences[0].index_match << endl;
  cout << "Correspondence 0 distance: " << all_correspondences[0].distance << endl;
  cout << "Correspondence 0 weight: " << all_correspondences[0].weight << endl;

  //------------------------
  // Determine the sign
  // init
  Eigen::Vector3f matchPoint;
  Eigen::Vector3f matchNormal;
  Eigen::Vector3f delta;
  Eigen::Vector3f queryVec;
  float dotProd;

  // Loop for each data point
  for (int j = 0; j < nDataPoints; j++){

    // Form the query vector
    queryVec(0) = query(0,all_correspondences[j].index_query);
    queryVec(1) = query(1,all_correspondences[j].index_query);
    queryVec(2) = query(2,all_correspondences[j].index_query);
    queryVec(0) = query_cloud->at(j).x;

    // cout << "QueryVec with index " << all_correspondences[j].index_query << " at j " << j << " is:\n" << queryVec << endl;

    // Get the point of the match and the normal at that point
    matchPoint(0) = cloud->at(all_correspondences[j].index_match).x;
    matchPoint(1) = cloud->at(all_correspondences[j].index_match).y;
    matchPoint(2) = cloud->at(all_correspondences[j].index_match).z;
    matchNormal(0) = -cloud->at(all_correspondences[j].index_match).normal_x;
    matchNormal(1) = -cloud->at(all_correspondences[j].index_match).normal_y;
    matchNormal(2) = -cloud->at(all_correspondences[j].index_match).normal_z;
    // TODO TDBM - work out why the normal needs to be negated

    // cout << "Normal is: " << matchNormal << endl;
    // cout << "Match point is:\n" << matchPoint << endl;

    // Get the vector from the surface to the test point 
    delta = queryVec - matchPoint;

    // cout << "Delta is:\n" << delta << endl;

    // Get the dot product - negative if the point is inside the surface (a dot b) = |a||b|cos(theta). |theta| > 90 deg means it is inside surface - results in a negative dot product
    dotProd = queryVec.dot(matchNormal);
    // cout << "Dot product result is: " << dotProd << endl;

    // Set the sign of the distance based on the dot product sign// TDBM buffer on signed check?
    if (dotProd < 0){
      // Point inside object
      distAndGrad(0,j) = std::sqrt(all_correspondences[j].distance)*-1.0;
      // cout << "Signed distance is: " << distances(j) << endl;
    }else{
      // Point Outside object
      distAndGrad(0,j) = std::sqrt(all_correspondences[j].distance);
      // cout << "Signed distance is: " << distances(j) << endl;
    }

    // Store the gradient (approxmated by using the vector from the matched surface point to the query point)
    distAndGrad(1,j) = delta(0);
    distAndGrad(2,j) = delta(1);
    distAndGrad(3,j) = delta(2);

  }

  return distAndGrad;
}
#include <ros/ros.h>
#include <iostream>

#include <nurbsS.h>
#include <cmath>

#include "cans/mapping3D.h"

#define PI 3.14159265

using namespace PLib;

Vector_HPoint3Df getRow(Matrix_HPoint3Df data,int row_id){
  
  int n_ctrl = data.rows();
  
  Vector_HPoint3Df ctrlRow(n_ctrl);
  
  for (int j = 0; j < n_ctrl; j++){
    ctrlRow[j] = data(row_id,j);
  }

  return ctrlRow;
}

Vector_HPoint3Df getCol(Matrix_HPoint3Df& data, int col_id){
  
  int n_ctrl = data.rows();
  
  Vector_HPoint3Df ctrlRow(n_ctrl);
  
  for (int i = 0; i < n_ctrl; i++){
    ctrlRow[i] = data(i,col_id);
  }

  return ctrlRow;
}

void insertRow(Matrix_HPoint3Df& data,Vector_HPoint3Df row, int row_id){
  // Function to insert data into an exisitng matrix - across a whole row
  for (int j = 0; j < data.cols(); j++){
    data(row_id,j) = row[j];
  }

}

void insertCol(Matrix_HPoint3Df& data,Vector_HPoint3Df col, int col_id){
  // Function to insert data into an exisitng matrix - across a whole row
  for (int i = 0; i < data.rows(); i++){
    data(i,col_id) = col[i];
  }

}


void testSurfaceJoin(){
  
  Mapping3D map;

  // Create surface
  int n_points = 11;
  Matrix_Point3Df scan1(n_points,n_points);
  Matrix_Point3Df scan2(n_points,n_points);
    
  for (int i = 0; i<scan1.rows(); i++){
    for (int j = 0; j<scan1.cols(); j++){
      scan1(i,j) = Point3Df((float)j/(n_points-1),(float)i/(n_points-1),0.0);

      scan2(i,j) = Point3Df((float)j/(n_points-1)+1.0,(float)i/(n_points-1),0.0);
    } 
  }
  
  cout << "scan (0,0) is: " << scan1(0,0) << endl;
  cout << "scan (end,0) is: " << scan1(n_points-1,0) << endl;
  cout << "scan (end,end) is: " << scan1(n_points-1,n_points-1) << endl;

  // Fit surfaces
  NurbsSurfacef srf1;
  NurbsSurfacef srf2;
  int deg = 3;
  int n_ctrl = 5;

  srf1.leastSquares(scan1,deg,deg,n_ctrl,n_ctrl);
  srf2.leastSquares(scan2,deg,deg,n_ctrl,n_ctrl);

  cout << "srf 1 (0,0) is: " << srf1(0,0) << endl;
  cout << "srf 1 (1,0) is: " << srf1(1.0,0.0) << endl;
  cout << "srf 1 (0,1) is: " << srf1(0.0,1.0) << endl;
  cout << "srf 1 (1,1) is: " << srf1(1.0,1.0) << endl;

  // Create a curve from a row of control points 
  Vector_HPoint3Df ctrlRow = getRow(srf1.ctrlPnts(),0);
  Vector_FLOAT knots = srf1.knotV();

  NurbsCurvef crvRow(ctrlRow,knots,deg);

  Vector_HPoint3Df ctrlRow2 = getRow(srf2.ctrlPnts(),0);
  Vector_FLOAT knots2 = srf2.knotV();

  NurbsCurvef crvRow2(ctrlRow2,knots2,deg);

  cout << "curve 1 ctrl pnts: " << crvRow.ctrlPnts() << endl;
  cout << "curve 2 ctrl pnts: " << crvRow2.ctrlPnts() << endl;

  // Join curves 
  NurbsCurvef crv_out = map.joinCurves(crvRow,crvRow2);

  cout << "curve out ctrl pnts: " << crv_out.ctrlPnts() << endl;

  Matrix_HPoint3Df ctrl_out(n_ctrl,2*n_ctrl - deg);
  Vector_FLOAT knots_out(2*n_ctrl + 1);
  Vector_HPoint3Df ctrlRow1;

  Vector_FLOAT knots1 = srf1.knotV();
  knots2 = srf2.knotV();

  for (int i = 0; i < n_ctrl; i++) {
    // Get rows out
    ctrlRow1 = getRow(srf1.ctrlPnts(),i);  
    ctrlRow2 = getRow(srf2.ctrlPnts(),i);
    

    // Create Curves
    NurbsCurvef crv1(ctrlRow1,knots1,deg);
    NurbsCurvef crv2(ctrlRow2,knots2,deg);

    // Join Curves
    crv_out = map.joinCurves(crv1,crv2);

    // Store surves in control point matrix
    insertRow(ctrl_out,crv_out.ctrlPnts(),i);


    knots_out += crv_out.knot();

  }

  // Average knot vector
  for (int i = 1; i < knots_out.size(); i++){
    knots_out[i] /= n_ctrl;
  }

  cout << "Output knots: " << knots_out << endl;
  cout << "output cntrl (0,0): " << ctrl_out(0,0) << endl;
  cout << "output cntrl (end,0): " << ctrl_out(n_ctrl-1,0) << endl;
  cout << "output cntrl (0,end): " << ctrl_out(0,2*n_ctrl - deg - 1) << endl;
  cout << "output cntrl (end,end): " << ctrl_out(n_ctrl-1,2*n_ctrl - deg - 1) << endl;

  // Create surface from new control points
  NurbsSurfacef srf3(deg,deg,srf1.knotU(),knots_out,ctrl_out);

  cout << "srf3 at (0,0): " << srf3(0.0,0.0) << endl;
  cout << "srf3 at (1.0,0): " << srf3(1.0,0.0) << endl;
  cout << "srf3 at (1.0,1.0): " << srf3(1.0,1.0) << endl;



  //int writePS(const char*, int nu, int nv, const Point_nD<T,N>& camera, const Point_nD<T,N>& lookAt, int cp=0,T magFact=T(-1),T dash=T(5)) const ;
  Point3Df cam(10,10,10);
  Point3Df look(0.0,0.0,0.0);
  srf1.writePS("srf1.ps",15,15,cam,look);
  srf2.writePS("srf2.ps",15,15,cam,look);
  srf3.writePS("srf3.ps",15,15,cam,look);

  srf3.writeVRML("tleastS.wrl",Color(255,100,255),50,80) ;
  // S.writePS("tleastS.ps",5,5,Point3Df(10,10,10),Point3Df(0,0,0)) ;

  char * extendDirection = "RU";

  cout << "extend Direction 1" << extendDirection[0] << endl;
  cout << "extend Direction 2" << extendDirection[1] << endl;

  if (extendDirection[0] == 'R'){// single quotations for the symbol
    cout << "Comparison with the if condition worked!";
  }

  cout << "strcmp check is : " << strcmp(extendDirection,"RU") << endl;


  NurbsSurfacef srf4 = map.joinSurfaces(srf1,srf2,"RR")

  srf4.writePS("srf4.ps",15,15,cam,look);
}



int main (int argc, char ** argv){

  ros::init (argc, argv, "surf_join");
  ros::NodeHandle nh;

  // FIt a NURBS surface
  testSurfaceJoin();

}
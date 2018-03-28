#include <ros/ros.h>
#include <iostream>

#include <nurbs.h>
#include <cmath>

#include "cans/mapping3D.h"

#define PI 3.14159265

using namespace PLib;

void testCurves(){

  // Setup measurements
  int n_param = 20;

  // Define the truth function 
  Vector_Point3Df obs1(20);
  Vector_Point3Df obs2(20);
  float t = 0.1;
  int j = 0;

  for (int i = 0; i< n_param*1.75; i++){
    t = 0.1 + i*(2*PI)/(n_param*1.75);
    if (i < 20){
      obs1[i].x() = sin(t);
      obs1[i].y() = pow(t,0.2) + cos(t);
      obs1[i].z() = 0.0;
    } 
    if (i >= 15){
      obs2[j].x() = sin(t) + 0.1;
      obs2[j].y() = pow(t,0.2) + cos(t) + 0.1;
      obs2[j].z() = 0.0;
      j++;
    }

  }
  
  cout << "Obs 1 are: " << obs1 << endl;

  // Create nurbs curve
  NurbsCurvef crv1;
  NurbsCurvef crv2;

  crv1.leastSquares(obs1,3,10);
  crv2.leastSquares(obs2,3,10);
  //int leastSquares(const Vector< Point_nD<T,N> >& Q, int degC, int n ) ;

  cout << "Curve1 control points are: " << crv1.ctrlPnts() << endl;
  cout << "Curve1 at (0.4) is: " << crv1.hpointAt(0.4) << endl;
  cout << "Curve1 at (0.4) is (with () operator): " << crv1(0.4) << endl;
  cout << "Curve1 4th knot is: " << crv1.knot(4) << endl;
  cout << "Curve1 degree is: " << crv1.degree() << endl;

  cout << "\n\nKnot Vector starts as: " << crv1.knot() << endl;
  cout << "Knot Vector 2 starts as: " << crv2.knot() << endl;

  // crv1.knot().push_back(1.1); // This doesn/t work as the push_back is not public

  Vector_FLOAT knots1(crv1.knot());

  // Trim curve 1 from 0:n_ctrl
  knots1.trim(crv1.ctrlPnts().size()+1);

  // Take 
  cout << "Size 1: " << crv2.knot().size() - crv2.degree() - 1 << endl;
  cout << "Size 2: " << crv2.ctrlPnts().size() << endl;
  // Vector_FLOAT knots2(crv2.knot().size() - crv2.degree());
  Vector_FLOAT knots2(crv2.ctrlPnts().size());

  for (int i = 0; i<crv2.ctrlPnts().size(); i++){
    knots2[i] = crv2.knot()[i+crv2.degree()+1] + crv1.maxKnot();
  }

  cout << "\n\nKnot Vector 1 step 2: " << knots1 << endl;
  cout << "Knot Vector 2 step 2: " << knots2 << endl;

  // Combine knot vectors together
  knots1.resize(2*crv1.ctrlPnts().size()+1);
  cout << "\n\nKnot Vector 1 step 3: " << knots1 << endl;

  j = 0;
  for (int i = 0; i < knots2.size(); i++ ){
    j = i + crv1.ctrlPnts().size()+1;

    knots1[j] = knots2[i];
  }

  cout << "\n\nKnot Vector 1 step 4: " << knots1 << endl;

  // Scale to 0 - 1
  for (int i = 0; i < knots1.size(); i++){
    knots1[i] /= knots1[knots1.size()-1];
  }
  // Would be more efficient to just divide what ever comes first by 2 right at the start when first looping through
  // This gets more difficult when they get flipped though...
  cout << "\n\nKnot Vector 1 step 5: " << knots1 << endl;

  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  // Now lets try to do this more efficiently for Control points 
  int n_ctrl = crv1.ctrlPnts().size();
  int n_ctrl_out = 2*n_ctrl-crv1.degree();
  Vector_HPoint3Df ctrl_out(n_ctrl_out);

  cout << "Number of control points: " << n_ctrl << endl;
  cout << "control points 2: " << crv2.ctrlPnts() << endl;
  cout << "control points 2 [1]: " << project(crv2.ctrlPnts()[1]) << endl;

  cout << "control out is: " << ctrl_out << endl;

  j = 0;
  for (int i = 0; i < n_ctrl_out; i++){
    if (i < n_ctrl){ 
      ctrl_out[i] = crv1.ctrlPnts()[i];
    }
    else {
      j = crv2.degree() + i - (n_ctrl);
      ctrl_out[i] = crv2.ctrlPnts()[j];
    }
    
  }

  cout << "control points combined " << ctrl_out << endl;

  // Create new NURBS curve
  // Vector_FLOAT weights(1.0,ctrl_out.size()); // Weights
  cout << "New ctl size: " << ctrl_out.size() << "\nNew Knots size: " << knots1.size() << endl;

  NurbsCurvef crv3(ctrl_out,knots1,crv1.degree());

  cout << "Curve 1 at param 0: " << crv1(0.0) << endl;
  cout << "Curve 3 at param 0: " << crv3(0.0) << endl;

  // Write to file to check
  crv1.writePS("crv1.ps");
  crv2.writePS("crv2.ps");
  crv3.writePS("crv3.ps");

  Mapping3D mp;

  NurbsCurvef crv4 = mp.joinCurves(crv1,crv2);
  cout << "curve 4 (A) knots: " << crv4.knot() << endl;

  NurbsCurvef crv5 = mp.joinCurves(crv1,crv2,false,true);

  crv4.writePS("crv4.ps");
  cout << "curve 5 (B) knots: " << crv5.knot() << endl;
  crv5.writePS("crv5.ps");

  NurbsCurvef crv6 = mp.joinCurves(crv2,crv1,true,false);
  cout << "curve 6 (C) knots: " << crv6.knot() << endl;

  NurbsCurvef crv7 = mp.joinCurves(crv1,crv2,true,true);
  cout << "curve 7 (C) knots: " << crv7.knot() << endl;


  crv6.writePS("crv6.ps");
  crv7.writePS("crv7.ps");
    
}


int main (int argc, char ** argv){
  // Initialize ROS
  ros::init (argc, argv, "test_mapping");
  ros::NodeHandle nh;

  // FIt a NURBS surface
  testCurves();

  // Spin
//   ros::spin ();
}

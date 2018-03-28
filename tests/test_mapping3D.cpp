#include <nurbs.h>
#include <cmath>

#include "cans/mapping3D.h"
#include "gtest/gtest.h"

#define PI 3.14159265

namespace {

using namespace PLib;

class mapping3DTest : public ::testing::Test {
 protected:
  mapping3DTest() : scan(3,3) {
    float k;
  
    for (int i = 0; i<scan.rows(); i++){
      for (int j=0; j<scan.cols(); j++){
        k = i*j + 2.0*i - j*3.0;
        scan(i,j) = Point3Df(i,j,k);
      } 
    }

  }

  ~mapping3DTest() {;}

  Matrix_Point3Df scan;


};

TEST_F (mapping3DTest, BlankConstructor){
    Mapping3D map;
}

TEST_F (mapping3DTest, TestComputeCentre){
    // Initialise class
    Mapping3D map;

    // Get centre
    Point3Df centre = map.compute_centre_of_data(scan);
    Point3Df expect(1,1,0);

    EXPECT_EQ(centre, expect);
}


//--------------------------------------------
//-------- Join Curve Tests -----
//--------------------------------------------

class nurbsJoiningTest : public ::testing::Test {
 protected:
  nurbsJoiningTest() : obs1(20), obs2(20), deg(3), n_ctrl(7) {
    
    int n_param = 20;
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

    crv1.leastSquares(obs1,deg,n_ctrl);
    crv2.leastSquares(obs1,deg,n_ctrl);
  }

  ~nurbsJoiningTest() {;}

  Vector_Point3Df obs1;
  Vector_Point3Df obs2;

  NurbsCurvef crv1;
  NurbsCurvef crv2;

  int deg;
  int n_ctrl;

};

TEST_F (nurbsJoiningTest, TestJoinCrv1Crv2){
    Mapping3D map;

    // Curve 1 before curve 2
    NurbsCurvef crv3 = map.joinCurves(crv1,crv2,false,false);

    // Check number of control points (n_ctrl1 + n_ctrl2 - degree)
    int n_ctrl_out = crv3.ctrlPnts().size();
    EXPECT_EQ(2*n_ctrl - deg,n_ctrl_out);

    // Check number of knots (n_ctrl1 + n_ctrl2 + 1) = n_ctrl_out + deg + 1
    int n_knot_out = crv3.knot().size();
    EXPECT_EQ(2*n_ctrl + 1,n_knot_out);

    // Check start and end knots - deg + 1 0's 
    for (int i = 0; i < deg+1; i ++){
        // Start knots
        EXPECT_FLOAT_EQ(0.0,crv3.knot()[i]);

        // end knots
        EXPECT_FLOAT_EQ(1.0,crv3.knot()[n_knot_out - 1 - i]);
    }

    // Check first control point 
    EXPECT_EQ(crv1.ctrlPnts()[0],crv3.ctrlPnts()[0]);

    // Check last control point 
    EXPECT_EQ(crv2.ctrlPnts()[n_ctrl-1],crv3.ctrlPnts()[n_ctrl_out-1]);

}

TEST_F (nurbsJoiningTest, TestJoinCrv1REVCrv2){
    Mapping3D map;

    // Curve 1 before curve 2, with curve 2 paramaters flipped
    NurbsCurvef crv3 = map.joinCurves(crv1,crv2,false,true);

    // Check number of control points (n_ctrl1 + n_ctrl2 - degree)
    int n_ctrl_out = crv3.ctrlPnts().size();
    EXPECT_EQ(2*n_ctrl - deg,n_ctrl_out);

    // Check number of knots (n_ctrl1 + n_ctrl2 + 1) = n_ctrl_out + deg + 1
    int n_knot_out = crv3.knot().size();
    EXPECT_EQ(2*n_ctrl + 1,n_knot_out);

    // Check start and end knots - deg + 1 0's 
    for (int i = 0; i < deg+1; i ++){
        // Start knots
        EXPECT_FLOAT_EQ(0.0,crv3.knot()[i]);

        // end knots
        EXPECT_FLOAT_EQ(1.0,crv3.knot()[n_knot_out - 1 - i]);
    }

    // Check first control point 
    EXPECT_EQ(crv1.ctrlPnts()[0],crv3.ctrlPnts()[0]);

    // Check last control point 
    EXPECT_EQ(crv2.ctrlPnts()[0],crv3.ctrlPnts()[n_ctrl_out-1]);

}

TEST_F (nurbsJoiningTest, TestJoinCrv2Crv1){
    Mapping3D map;

    // Curve 2 before curve 1
    NurbsCurvef crv3 = map.joinCurves(crv1,crv2,true,false);

    // Check number of control points (n_ctrl1 + n_ctrl2 - degree)
    int n_ctrl_out = crv3.ctrlPnts().size();
    EXPECT_EQ(2*n_ctrl - deg,n_ctrl_out);

    // Check number of knots (n_ctrl1 + n_ctrl2 + 1) = n_ctrl_out + deg + 1
    int n_knot_out = crv3.knot().size();
    EXPECT_EQ(2*n_ctrl + 1,n_knot_out);

    // Check start and end knots - deg + 1 0's 
    for (int i = 0; i < deg+1; i ++){
        // Start knots
        EXPECT_FLOAT_EQ(0.0,crv3.knot()[i]);

        // end knots
        EXPECT_FLOAT_EQ(1.0,crv3.knot()[n_knot_out - 1 - i]);
    }

    // Check first control point 
    EXPECT_EQ(crv2.ctrlPnts()[0],crv3.ctrlPnts()[0]);

    // Check last control point 
    EXPECT_EQ(crv1.ctrlPnts()[n_ctrl-1],crv3.ctrlPnts()[n_ctrl_out-1]);

}

TEST_F (nurbsJoiningTest, TestJoinREVCrv2Crv1){
    Mapping3D map;

    // Curve 2 before curve 1, with curve 2 paramaters flipped
    NurbsCurvef crv3 = map.joinCurves(crv1,crv2,true,true);

    // Check number of control points (n_ctrl1 + n_ctrl2 - degree)
    int n_ctrl_out = crv3.ctrlPnts().size();
    EXPECT_EQ(2*n_ctrl - deg,n_ctrl_out);

    // Check number of knots (n_ctrl1 + n_ctrl2 + 1) = n_ctrl_out + deg + 1
    int n_knot_out = crv3.knot().size();
    EXPECT_EQ(2*n_ctrl + 1,n_knot_out);

    // Check start and end knots - deg + 1 0's 
    for (int i = 0; i < deg+1; i ++){
        // Start knots
        EXPECT_FLOAT_EQ(0.0,crv3.knot()[i]);

        // end knots
        EXPECT_FLOAT_EQ(1.0,crv3.knot()[n_knot_out - 1 - i]);
    }

    // Check first control point 
    EXPECT_EQ(crv2.ctrlPnts()[n_ctrl-1],crv3.ctrlPnts()[0]);

    // Check last control point 
    EXPECT_EQ(crv1.ctrlPnts()[n_ctrl-1],crv3.ctrlPnts()[n_ctrl_out-1]);

}

TEST_F (nurbsJoiningTest, TestJoinDifferentNctrl){
    Mapping3D map;

    // Change the number of control points in curve 2
    int n_ctrl2 = 10;

    crv2.leastSquares(obs2,deg,n_ctrl2);

    // Curve 2 before curve 1, with curve 2 paramaters flipped
    NurbsCurvef crv3 = map.joinCurves(crv1,crv2,false,false);

    // Check number of control points (n_ctrl1 + n_ctrl2 - degree)
    int n_ctrl_out = crv3.ctrlPnts().size();
    EXPECT_EQ(n_ctrl + n_ctrl2 - deg,n_ctrl_out);

    // Check number of knots (n_ctrl1 + n_ctrl2 + 1) = n_ctrl_out + deg + 1
    int n_knot_out = crv3.knot().size();
    EXPECT_EQ(n_ctrl + n_ctrl2 + 1,n_knot_out);

    // Check start and end knots - deg + 1 0's 
    for (int i = 0; i < deg+1; i ++){
        // Start knots
        EXPECT_FLOAT_EQ(0.0,crv3.knot()[i]);

        // end knots
        EXPECT_FLOAT_EQ(1.0,crv3.knot()[n_knot_out - 1 - i]);
    }

    // Check first control point 
    EXPECT_EQ(crv1.ctrlPnts()[0],crv3.ctrlPnts()[0]);

    // Check last control point 
    EXPECT_EQ(crv2.ctrlPnts()[n_ctrl2-1],crv3.ctrlPnts()[n_ctrl_out-1]);

}

} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
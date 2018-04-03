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

//--------------------------------------------
//-------- Join Surface Tests -----
//--------------------------------------------

class nurbsSurfaceJoiningTest : public ::testing::Test {
 protected:
  nurbsSurfaceJoiningTest() : deg(3), n_ctrl(5) {
    
    int n_points = 11;
    
    Matrix_Point3Df scan1A(n_points,n_points);
    Matrix_Point3Df scan1B(n_points,n_points);
    Matrix_Point3Df scan1C(n_points,n_points);
    Matrix_Point3Df scan1D(n_points,n_points);        
    Matrix_Point3Df scan2A(n_points,n_points);
    Matrix_Point3Df scan2B(n_points,n_points);
    Matrix_Point3Df scan2C(n_points,n_points);
    Matrix_Point3Df scan2D(n_points,n_points);

    for (int i = 0; i<n_points; i++){
        for (int j = 0; j<n_points; j++){
            scan1A(i,j) = Point3Df((float)j/(n_points-1),(float)i/(n_points-1),0.0);
            scan1B(i,j) = Point3Df((float)j/(n_points-1) + 2.0,(float)i/(n_points-1),0.0);
            scan1C(i,j) = Point3Df((float)j/(n_points-1) + 1.0,(float)i/(n_points-1) - 1.0,0.0);
            scan1D(i,j) = Point3Df((float)j/(n_points-1) + 1.0,(float)i/(n_points-1) + 1.0,0.0);
            
            
            scan2A(i,j) = Point3Df((float)j/(n_points-1)+1.0,(float)i/(n_points-1),0.0);
            scan2B(n_points - 1 - i,n_points - 1 - j) = Point3Df((float)j/(n_points-1)+1.0,(float)i/(n_points-1),0.0);
            scan2C(j,n_points - 1 - i) = Point3Df((float)j/(n_points-1)+1.0,(float)i/(n_points-1),0.0);
            scan2D(n_points - 1 - j,i) = Point3Df((float)j/(n_points-1)+1.0,(float)i/(n_points-1),0.0);
            // cout << "(i, j) = (" << n_points - 1 - j << ", " << i << "), scan2D is: " << scan2D(n_points - 1 - j,i) << endl;
        } 

        
    }
    // cout << "Scan2C (1,1): " << scan2C(10,10) << endl;
    // for (int i = 0; i < 4; i ++){
    //     srf1_arr[i]* = new NurbsSurfacef;
    //     srf2_arr[i]* = new NurbsSurfacef;

    //     srf1_arr[i].leastSquares(scan1A,deg,deg,n_ctrl,n_ctrl);
    //     srf2_arr[i].leastSquares(scan2A,deg,deg,n_ctrl,n_ctrl);
    // }

    srf1A.leastSquares(scan1A, deg, deg, n_ctrl, n_ctrl);
    srf1B.leastSquares(scan1B, deg, deg, n_ctrl, n_ctrl);
    srf1C.leastSquares(scan1C, deg, deg, n_ctrl, n_ctrl);
    srf1D.leastSquares(scan1D, deg, deg, n_ctrl, n_ctrl);

    srf2A.leastSquares(scan2A, deg, deg, n_ctrl, n_ctrl);
    srf2B.leastSquares(scan2B, deg, deg, n_ctrl, n_ctrl);
    srf2C.leastSquares(scan2C, deg, deg, n_ctrl, n_ctrl);
    srf2D.leastSquares(scan2D, deg, deg, n_ctrl, n_ctrl);
    
  }

  
//   NurbsSurfaceArray srf1_arr(NurbsSurfacef * srf1, 4);
//   NurbsSurfaceArray srf2_arr(NurbsSurfacef * srf2, 4);

  NurbsSurfacef srf1A;
  NurbsSurfacef srf1B;
  NurbsSurfacef srf1C;
  NurbsSurfacef srf1D;
  NurbsSurfacef srf2A;
  NurbsSurfacef srf2B;
  NurbsSurfacef srf2C;
  NurbsSurfacef srf2D;

  int deg;
  int n_ctrl;

  Mapping3D map;

};


TEST_F (nurbsSurfaceJoiningTest, testRR){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1A,srf2A,"RR");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(0.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(0.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}

TEST_F (nurbsSurfaceJoiningTest, testRL){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1A,srf2B,"RL");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(0.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(0.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}

TEST_F (nurbsSurfaceJoiningTest, testRU){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1A,srf2C,"RU");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(0.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(0.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}

TEST_F (nurbsSurfaceJoiningTest, testRD){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1A,srf2D,"RD");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(0.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(0.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}

//------------------------------------------
TEST_F (nurbsSurfaceJoiningTest, testLL){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1B,srf2A,"LL");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(3.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(3.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}

TEST_F (nurbsSurfaceJoiningTest, testLR){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1B,srf2B,"LR");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(3.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(3.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}

TEST_F (nurbsSurfaceJoiningTest, testLD){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1B,srf2C,"LD");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(3.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(3.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}

TEST_F (nurbsSurfaceJoiningTest, testLU){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1B,srf2D,"LU");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotV().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(3.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(3.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());
}
// --------------------------------------
// Move srf 1 to be below srf 2
// --------------------------------------
TEST_F (nurbsSurfaceJoiningTest, testUU){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1C,srf2A,"UU");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,-1.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,-1.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}

TEST_F (nurbsSurfaceJoiningTest, testUD){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1C,srf2B,"UD");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,-1.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,-1.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}

TEST_F (nurbsSurfaceJoiningTest, testUL){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1C,srf2C,"UL");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,-1.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,-1.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}

TEST_F (nurbsSurfaceJoiningTest, testUR){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1C,srf2D,"UR");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,-1.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,-1.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,1.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}
// --------------------------------------
// Move srf 1 to be above srf 2
// --------------------------------------
TEST_F (nurbsSurfaceJoiningTest, testDD){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1D,srf2A,"DD");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,2.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,2.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}

TEST_F (nurbsSurfaceJoiningTest, testDU){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1D,srf2B,"DU");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,2.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,2.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}

TEST_F (nurbsSurfaceJoiningTest, testDR){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1D,srf2C,"DR");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,2.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,2.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}

TEST_F (nurbsSurfaceJoiningTest, testDL){

    NurbsSurfacef srfOut = map.joinSurfaces(srf1D,srf2D,"DL");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1,srfOut.knotU().size());

    // Test Corner points 
    EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),srfOut(0.0,0.0));
    EXPECT_EQ(HPoint3Df(1.0,2.0,0.0,1.0),srfOut(1.0,0.0));
    EXPECT_EQ(HPoint3Df(2.0,0.0,0.0,1.0),srfOut(0.0,1.0));
    EXPECT_EQ(HPoint3Df(2.0,2.0,0.0,1.0),srfOut(1.0,1.0));

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());   
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());
}


class meshFromScanTest : public ::testing::Test{
    protected:
        meshFromScanTest() : cloud(new pcl::PointCloud<pcl::PointXYZ>(10,10)), nanArray(10,10), expectNanArray(10,10),
                rowFlags(10,1), colFlags(1,10)
        {
        
            int n_points = 10;
            bool nanFlag = false;

            // Set to all zeros
            expectNanArray.setZero();
            nanArray.setZero();

            // Create data with Nans
            for (int i = 0; i<n_points; i++){
                for (int j = 0; j<n_points; j++){
                    if (i <= 1 || (i == 2 && (j < 3 || j > 7)) || i == 9 || j == 0 || j == 9 || (i == 5 && j == 4)){
                        nanFlag = true;
                    }
                    if (nanFlag) {
                        cloud->at(j,i).x = NAN;
                        cloud->at(j,i).y = NAN;
                        cloud->at(j,i).z = NAN;
                        expectNanArray(i,j) = true;
                    }else{
                        cloud->at(j,i).x = (float)j/(n_points-1);
                        cloud->at(j,i).y = (float)i/(n_points-1);
                        cloud->at(j,i).z = 0.0;//(float)i*2.0f - (float)j*5.0f;
                    }
                    nanFlag = false;
                } 
            }

            // Set flags to 1s
            rowFlags.setOnes(10, 1);
            colFlags.setOnes(1, 10);

            // Set to remove all nans
            mp.maxNanAllowed = 0;
            mp.numRowsDesired = 2;
            mp.numColsDesired = 2;

        }

        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud;

        Eigen::Array<bool, 10, 10> expectNanArray;

        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> nanArray;

        // Initialise
        Eigen::Array<bool, Eigen::Dynamic, 1> rowFlags;
        Eigen::Array<bool, 1, Eigen::Dynamic> colFlags;
  
        // Init class
        Mapping3D mp;
};

TEST_F (meshFromScanTest, nanMatrixFromPCTest){

  // Fill the Nan Array
  mp.getNanMatrixFromPointCloud(nanArray, cloud);

  cout << "Expected Nan array is:\n" << expectNanArray << endl;
  cout << "Computed Nan array is:\n" << nanArray << endl;

  EXPECT_TRUE((nanArray==expectNanArray).all());

}

TEST_F (meshFromScanTest, rowNanRemoveOnce){


    // Fill the Nan Array
    mp.getNanMatrixFromPointCloud(nanArray, cloud);
    
    // Run the function
    int nansPresent = mp.removeRowsNan(nanArray, rowFlags);
    
    EXPECT_TRUE(nansPresent);

    Eigen::Array<bool, 10, 1> expectRowFlags(10,1); expectRowFlags.setOnes(10,1);
    expectRowFlags[0] = false;
    expectRowFlags[1] = false;
    expectRowFlags[9] = false;

    cout << "Expected Row Flags: " << expectRowFlags << endl;
    cout << "Computed Row Flags: " << rowFlags << endl;

    EXPECT_TRUE((expectRowFlags==rowFlags).all());

    EXPECT_EQ(false,nanArray(0,0));
    EXPECT_EQ(false,nanArray(1,1));
    EXPECT_EQ(false,nanArray(9,9));
}

TEST_F (meshFromScanTest, colNanRemoveOnce){


    // Fill the Nan Array
    mp.getNanMatrixFromPointCloud(nanArray, cloud);
    
    // Run the function
    int nansPresent = mp.removeColsNan(nanArray, colFlags);
    
    EXPECT_TRUE(nansPresent);

    Eigen::Array<bool, 1, 10> expectColFlags(1,10); expectColFlags.setOnes(1,10);
    expectColFlags[0] = false;
    expectColFlags[9] = false;

    cout << "Expected Col Flags: " << expectColFlags << endl;
    cout << "Computed Col Flags: " << colFlags << endl;

    EXPECT_TRUE((expectColFlags==colFlags).all());
    EXPECT_EQ(false,nanArray(0,0));
    EXPECT_EQ(false,nanArray(0,9));
}

TEST_F (meshFromScanTest, colRowNanRemove){


    // Fill the Nan Array
    mp.getNanMatrixFromPointCloud(nanArray, cloud);
    
    bool exitFlag = false;

    bool nansPresent = true;

    // Loop to remove Nans
    while (!exitFlag){
        nansPresent = mp.removeRowsNan(nanArray, rowFlags);

        if (nansPresent){
        nansPresent = mp.removeColsNan(nanArray, colFlags);
        if (nansPresent){
            // Check dimensions 
            cout << "Nans removed: Row count: " << rowFlags.count() << "\nCol Count: " << colFlags.count() << endl;
            if (rowFlags.count() <= 0 || colFlags.count() <= 0){
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
    
    Eigen::Array<bool, 10, 1> expectRowFlags(10,1); expectRowFlags.setOnes(10,1);
    expectRowFlags[0] = false;
    expectRowFlags[1] = false;
    expectRowFlags[2] = false;
    expectRowFlags[9] = false;

    Eigen::Array<bool, 1, 10> expectColFlags(1,10); expectColFlags.setOnes(1,10);
    expectColFlags[0] = false;
    expectColFlags[4] = false;
    expectColFlags[9] = false;

    cout << "Expected Row Flags: " << expectRowFlags << endl;
    cout << "Computed Row Flags: " << rowFlags << endl;

    EXPECT_TRUE((expectRowFlags==rowFlags).all());

    cout << "Expected Col Flags: " << expectColFlags << endl;
    cout << "Computed Col Flags: " << colFlags << endl;

    EXPECT_TRUE((expectColFlags==colFlags).all());

    EXPECT_EQ(false,nanArray(5,4));
}

TEST_F (meshFromScanTest, downsampleRowTestAllOnes){
    
    Eigen::Array<bool, Eigen::Dynamic, 1> rowFlagsT(10,1);rowFlagsT.setOnes(10,1);

    mp.numRowsDesired = 5;

    mp.downsampleRow(rowFlagsT);

    EXPECT_EQ(mp.numRowsDesired,rowFlagsT.count());

    cout << "rowFlags out: " << rowFlagsT << endl;

    EXPECT_EQ(true,rowFlagsT[0]);
    EXPECT_EQ(true,rowFlagsT[9]);
    EXPECT_EQ(false,rowFlagsT[1]);
    EXPECT_EQ(false,rowFlagsT[8]);
}

TEST_F (meshFromScanTest, downsampleRowTestMixed1){
    
    Eigen::Array<bool, Eigen::Dynamic, 1> rowFlagsT(10,1);rowFlagsT.setOnes(10,1);
    rowFlagsT[0] = false;
    rowFlagsT[1] = false;
    rowFlagsT[5] = false;
    rowFlagsT[6] = false;
    rowFlagsT[8] = false;

    mp.numRowsDesired = 5;

    mp.downsampleRow(rowFlagsT);

    EXPECT_EQ(mp.numRowsDesired,rowFlagsT.count());

    cout << "rowFlags out: " << rowFlagsT << endl;

    EXPECT_EQ(true,rowFlagsT[2]);
    EXPECT_EQ(true,rowFlagsT[9]);
    EXPECT_EQ(false,rowFlagsT[1]);
    EXPECT_EQ(false,rowFlagsT[8]);
}

TEST_F (meshFromScanTest, downsampleRowTestMixed2){
    
    Eigen::Array<bool, Eigen::Dynamic, 1> rowFlagsT(10,1);rowFlagsT.setOnes(10,1);
    rowFlagsT[0] = false;
    rowFlagsT[1] = false;
    rowFlagsT[5] = false;
 
    mp.numRowsDesired = 5;

    mp.downsampleRow(rowFlagsT);

    EXPECT_EQ(mp.numRowsDesired,rowFlagsT.count());

    cout << "rowFlags out: " << rowFlagsT << endl;

    EXPECT_EQ(true,rowFlagsT[2]);
    EXPECT_EQ(true,rowFlagsT[9]);
    EXPECT_EQ(false,rowFlagsT[1]);
}

TEST_F (meshFromScanTest, downsampleColTestAllOnes){
    
    Eigen::Array<bool, 1, Eigen::Dynamic> colFlagsT(1,10);colFlagsT.setOnes(1,10);

    mp.numColsDesired = 5;

    mp.downsampleCol(colFlagsT);

    EXPECT_EQ(mp.numColsDesired,colFlagsT.count());

    cout << "colFlags out: " << colFlagsT << endl;

    EXPECT_EQ(true,colFlagsT[0]);
    EXPECT_EQ(true,colFlagsT[9]);
    EXPECT_EQ(false,colFlagsT[1]);
    EXPECT_EQ(false,colFlagsT[8]);
}

TEST_F (meshFromScanTest, downsampleColTestMixed1){
    
    Eigen::Array<bool, 1, Eigen::Dynamic> colFlagsT(1,10);colFlagsT.setOnes(1,10);
    colFlagsT[0] = false;
    colFlagsT[1] = false;
    colFlagsT[5] = false;
    colFlagsT[6] = false;
    colFlagsT[8] = false;

    mp.numColsDesired = 5;

    mp.downsampleCol(colFlagsT);

    EXPECT_EQ(mp.numColsDesired,colFlagsT.count());

    cout << "colFlags out: " << colFlagsT << endl;

    EXPECT_EQ(true,colFlagsT[2]);
    EXPECT_EQ(true,colFlagsT[9]);
    EXPECT_EQ(false,colFlagsT[1]);
    EXPECT_EQ(false,colFlagsT[8]);
}

TEST_F (meshFromScanTest, downsampleColTestMixed2){
    
    Eigen::Array<bool, 1, Eigen::Dynamic> colFlagsT(1,10);colFlagsT.setOnes(1,10);
    colFlagsT[0] = false;
    colFlagsT[1] = false;
    colFlagsT[5] = false;
 
    mp.numColsDesired = 5;

    mp.downsampleCol(colFlagsT);

    EXPECT_EQ(mp.numColsDesired,colFlagsT.count());

    cout << "colFlags out: " << colFlagsT << endl;

    EXPECT_EQ(true,colFlagsT[2]);
    EXPECT_EQ(true,colFlagsT[9]);
    EXPECT_EQ(false,colFlagsT[1]);
}

// TEST_F (meshFromScanTest, regionAverageTest){

// }

TEST_F (meshFromScanTest, averageOutNanTest){

    Eigen::Array<int, 2, Eigen::Dynamic> nanIndices(2,1);
    nanIndices(0,0) = 5;
    nanIndices(1,0) = 4;

    // Average Nans
    mp.averageOutNans(cloud, nanIndices);

    pcl::PointXYZ expectP;

    expectP.x = (float)nanIndices(1,0)/(10-1);
    expectP.y = (float)nanIndices(0,0)/(10-1);
    expectP.z = 0.0;

    EXPECT_NEAR(expectP.x,cloud->at(nanIndices(1,0),nanIndices(0,0)).x,0.0001);
    EXPECT_NEAR(expectP.y,cloud->at(nanIndices(1,0),nanIndices(0,0)).y,0.0001);
    EXPECT_NEAR(expectP.z,cloud->at(nanIndices(1,0),nanIndices(0,0)).z,0.0001);

    cout << "Expected point:\n" << expectP << endl;
    cout << "Computed point:\n" << cloud->at(nanIndices(1,0),nanIndices(0,0)) << endl;
}

TEST_F (meshFromScanTest, meshFromScanFunction){

    mp.numRowsDesired = 6;
    mp.numColsDesired = 6;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOut (new pcl::PointCloud<pcl::PointXYZ>(mp.numColsDesired,mp.numRowsDesired,pcl::PointXYZ(0.0,0.0,0.0)));

    // Mesh from scan
    mp.meshFromScan(cloudOut, cloud);

    EXPECT_EQ(6,cloudOut->width);
    EXPECT_EQ(6,cloudOut->height);

    // Test the desired corner points
    EXPECT_FLOAT_EQ(cloud->at(1,3).x,cloudOut->at(0,0).x);
    EXPECT_FLOAT_EQ(cloud->at(1,3).y,cloudOut->at(0,0).y);
    EXPECT_FLOAT_EQ(cloud->at(1,3).z,cloudOut->at(0,0).z);

    EXPECT_FLOAT_EQ(cloud->at(1,8).x,cloudOut->at(0,5).x);
    EXPECT_FLOAT_EQ(cloud->at(1,8).y,cloudOut->at(0,5).y);
    EXPECT_FLOAT_EQ(cloud->at(1,8).z,cloudOut->at(0,5).z);
    
    EXPECT_FLOAT_EQ(cloud->at(8,8).x,cloudOut->at(5,5).x);
    EXPECT_FLOAT_EQ(cloud->at(8,8).y,cloudOut->at(5,5).y);
    EXPECT_FLOAT_EQ(cloud->at(8,8).z,cloudOut->at(5,5).z);

    EXPECT_FLOAT_EQ(cloud->at(8,3).x,cloudOut->at(5,0).x);
    EXPECT_FLOAT_EQ(cloud->at(8,3).y,cloudOut->at(5,0).y);
    EXPECT_FLOAT_EQ(cloud->at(8,3).z,cloudOut->at(5,0).z);

}

class meshFromScanLargeTest : public ::testing::Test{
    protected:
        meshFromScanLargeTest() : cloud(new pcl::PointCloud<pcl::PointXYZ>(1000,1000))
        {
        
            int n_points = 1000;
            bool nanFlag = false;

            // Create data with Nans
            for (int i = 0; i<n_points; i++){
                for (int j = 0; j<n_points; j++){
                    if (i <= 1 || (i == 2 && (j < 3 || j > 7)) || i == 9 || j == 0 || i == 5 || j == 578 || i == 345 || j == 9 || (i == 5 && j == 4)){
                        nanFlag = true;
                    }
                    if ()
                    if (nanFlag) {
                        cloud->at(j,i).x = NAN;
                        cloud->at(j,i).y = NAN;
                        cloud->at(j,i).z = NAN;
                    }else{
                        cloud->at(j,i).x = (float)j/(n_points-1);
                        cloud->at(j,i).y = (float)i/(n_points-1);
                        cloud->at(j,i).z = 0.0;//(float)i*2.0f - (float)j*5.0f;
                    }
                    nanFlag = false;
                } 
            }

            // Set to remove all nans
            mp.maxNanAllowed = 10;
            mp.numRowsDesired = 100;
            mp.numColsDesired = 100;

        }

        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud;
  
        // Init class
        Mapping3D mp;
};

TEST_F (meshFromScanLargeTest, meshFromScanFunctionLarge){

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOut (new pcl::PointCloud<pcl::PointXYZ>(mp.numColsDesired,mp.numRowsDesired,pcl::PointXYZ(0.0,0.0,0.0)));

    // Mesh from scan
    mp.meshFromScan(cloudOut, cloud);

    EXPECT_EQ(100,cloudOut->width);
    EXPECT_EQ(100,cloudOut->height);
}

} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
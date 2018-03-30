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

} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
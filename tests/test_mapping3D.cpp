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
//-------- NURBS Fitting tests  -----
//--------------------------------------------
class nurbsFittingTest : public ::testing::Test {
 protected:
  nurbsFittingTest() : obs1(20), obs2(20), obsSurf(20,20), deg(3), n_ctrl(7) {
    
    int n_param = 20;
    float t = 0.1;
    int j = 0;

    Eigen::Array<float,1,Eigen::Dynamic> vals(1,n_param);
    vals << 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,3,4;
    vals << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,2,3,3,4;

    for (int i = 0; i< n_param; i++){
        obs1[i].x() = (float)i/(float)n_param;
        obs1[i].y() = 0.0;
        obs1[i].z() = 0.0;

        // Multiple double-ups
        obs2[i].x() = vals(i);
        obs2[i].y() = 0.0;
        obs2[i].z() = 0.0;
        for (int j = 0; j < n_param; j++){
            obsSurf(i,j).x() = (float)j/(float)n_param;
            obsSurf(i,j).y() = vals(i);
            obsSurf(i,j).z() = 0.0;
        }
    }

    // Double up one
    obs1[0].x() = obs1[1].x();
    // crv1.leastSquares(obs1,deg,n_ctrl);
    // crv2.leastSquares(obs1,deg,n_ctrl);
  }

  ~nurbsFittingTest() {;}

  Vector_Point3Df obs1;
  Vector_Point3Df obs2;

  Matrix_Point3Df obsSurf;

  NurbsCurvef crv1;
  NurbsCurvef crv2;
  Object3D obj;

  int deg;
  int n_ctrl;

};

TEST_F (nurbsFittingTest, oneDoubleUpPoint){

    // Test that this doesn't give nans
    crv1.leastSquares(obs1,deg,n_ctrl);
    bool bAreThereNans = false;
    for (float s = 0.0; s < 1.0; s = s + 0.001){
        if (!std::isfinite(crv1(s).x())){
            bAreThereNans = true;
        }
    }

    EXPECT_TRUE(!bAreThereNans);
}

TEST_F (nurbsFittingTest, multipleDoubleUpPoints){

    // Test that this doesn't give nans
    crv2.leastSquares(obs2,deg,n_ctrl);
    bool bAreThereNans = false;
    for (float s = 0.0; s < 1.0; s = s + 0.001){
        if (!std::isfinite(crv2(s).x())){
            bAreThereNans = true;
        }
        if (!std::isfinite(crv2(s).y())){
            bAreThereNans = true;
        }
        if (!std::isfinite(crv2(s).z())){
            bAreThereNans = true;
        }
    }

    cout << "Curve Knots:\n" << crv2.knot() << endl;

    cout << "Control Points:\n" << crv2.ctrlPnts() << endl;

    // THIS SHOULD FAIL _ JUST CHECKING IT DOES
    EXPECT_TRUE(bAreThereNans);
}

TEST_F (nurbsFittingTest, nCtrlMoreNPoints){

    // Test that this doesn't give nans
    crv1.leastSquares(obs1,deg,21);
    bool bAreThereNans = false;
    for (float s = 0.0; s < 1.0; s = s + 0.001){
        if (!std::isfinite(crv1(s).x())){
            bAreThereNans = true;
        }
    }

    EXPECT_TRUE(!bAreThereNans);
}
/*
TEST_F (nurbsFittingTest, multipleDoubleUpPointsSurf){

    // Test that this doesn't give nans
    obj.leastSquares(obsSurf,deg,deg,n_ctrl,n_ctrl);
    bool bAreThereNans = false;
    for (float s = 0.0; s < 1.0; s = s + 0.01){
        for (float t = 0.0; t < 1.0; t = t + 0.01){
            if (!std::isfinite(obj(s,t).x())){
                bAreThereNans = true;
                cout << "NANS!" << endl;
            }
        }
    }

    EXPECT_TRUE(!bAreThereNans);
}

TEST_F (nurbsFittingTest, multipleDoubleUpPointsSurfConstructor){

    // Test that this doesn't give nans
    Object3D obj2(obsSurf,deg,deg,n_ctrl,n_ctrl);
    bool bAreThereNans = false;
    for (float s = 0.0; s < 1.0; s = s + 0.01){
        for (float t = 0.0; t < 1.0; t = t + 0.01){
            if (!std::isfinite(obj2(s,t).x())){
                bAreThereNans = true;
                cout << "NANS!" << endl;
            }
        }
    }

    EXPECT_TRUE(!bAreThereNans);
}
*/



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
/*
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
*/
//--------------------------------------------
//-------- Join Object3D Tests -----
//--------------------------------------------

class object3DJoiningTest : public ::testing::Test {
 protected:
  object3DJoiningTest() : deg(3), n_ctrl(5) {
    
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
    //     srf1_arr[i]* = new Object3D;
    //     srf2_arr[i]* = new Object3D;

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

  
//   NurbsSurfaceArray srf1_arr(Object3D * srf1, 4);
//   NurbsSurfaceArray srf2_arr(Object3D * srf2, 4);

  Object3D srf1A;
  Object3D srf1B;
  Object3D srf1C;
  Object3D srf1D;
  Object3D srf2A;
  Object3D srf2B;
  Object3D srf2C;
  Object3D srf2D;

  int deg;
  int n_ctrl;

  Mapping3D map;

};


TEST_F (object3DJoiningTest, testRR){

    Object3D srfOut = map.joinSurfaces(srf1A,srf2A,"RR");

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

    // srfOut.writeVRML("RR.wrl",Color(255,100,255),50,80);
}

TEST_F (object3DJoiningTest, testRL){

    Object3D srfOut = map.joinSurfaces(srf1A,srf2B,"RL");

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

TEST_F (object3DJoiningTest, testRU){

    Object3D srfOut = map.joinSurfaces(srf1A,srf2C,"RU");

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

TEST_F (object3DJoiningTest, testRD){

    Object3D srfOut = map.joinSurfaces(srf1A,srf2D,"RD");

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
TEST_F (object3DJoiningTest, testLL){

    Object3D srfOut = map.joinSurfaces(srf1B,srf2A,"LL");

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

TEST_F (object3DJoiningTest, testLR){

    Object3D srfOut = map.joinSurfaces(srf1B,srf2B,"LR");

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

TEST_F (object3DJoiningTest, testLD){

    Object3D srfOut = map.joinSurfaces(srf1B,srf2C,"LD");

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

    // srfOut.writeVRML("LD.wrl",Color(255,100,255),50,80);
}

TEST_F (object3DJoiningTest, testLU){

    Object3D srfOut = map.joinSurfaces(srf1B,srf2D,"LU");

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
TEST_F (object3DJoiningTest, testUU){

    Object3D srfOut = map.joinSurfaces(srf1C,srf2A,"UU");

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

TEST_F (object3DJoiningTest, testUD){

    Object3D srfOut = map.joinSurfaces(srf1C,srf2B,"UD");

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

TEST_F (object3DJoiningTest, testUL){

    Object3D srfOut = map.joinSurfaces(srf1C,srf2C,"UL");

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

    // srfOut.writeVRML("UL.wrl",Color(255,100,255),50,80);
}

TEST_F (object3DJoiningTest, testUR){

    Object3D srfOut = map.joinSurfaces(srf1C,srf2D,"UR");

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
TEST_F (object3DJoiningTest, testDD){

    Object3D srfOut = map.joinSurfaces(srf1D,srf2A,"DD");

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

TEST_F (object3DJoiningTest, testDU){

    Object3D srfOut = map.joinSurfaces(srf1D,srf2B,"DU");

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

    // srfOut.writeVRML("DU.wrl",Color(255,100,255),50,80);
}

TEST_F (object3DJoiningTest, testDR){

    Object3D srfOut = map.joinSurfaces(srf1D,srf2C,"DR");

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

TEST_F (object3DJoiningTest, testDL){

    Object3D srfOut = map.joinSurfaces(srf1D,srf2D,"DL");

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

TEST_F (object3DJoiningTest, testAsymRR){

    // Insert knots into the first surface
    Vector_FLOAT params(3);
    params[0] = 0.1;
    params[1] = 0.5;
    params[2] = 0.6;

    srf1A.refineKnotV(params);

    Object3D srfOut = map.joinSurfaces(srf1A,srf2A,"RR");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().rows());
    EXPECT_EQ(2*n_ctrl - deg + 3,srfOut.ctrlPnts().cols());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotU().size());
    EXPECT_EQ(2*n_ctrl + 1 + 3,srfOut.knotV().size());

    // Test Corner points 
    for (int i = 0; i < 4; i++){
        EXPECT_FLOAT_EQ(HPoint3Df(0.0,0.0,0.0,1.0).data[i],srfOut(0.0,0.0).data[i]);
        EXPECT_FLOAT_EQ(HPoint3Df(0.0,1.0,0.0,1.0).data[i],srfOut(1.0,0.0).data[i]);
        EXPECT_FLOAT_EQ(HPoint3Df(2.0,0.0,0.0,1.0).data[i],srfOut(0.0,1.0).data[i]);
        EXPECT_FLOAT_EQ(HPoint3Df(2.0,1.0,0.0,1.0).data[i],srfOut(1.0,1.0).data[i]);
    }
    
    // Test mid-extension points 
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.3).y());
    EXPECT_FLOAT_EQ(0.0,srfOut(0.0,0.7).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.3).y());
    EXPECT_FLOAT_EQ(1.0,srfOut(1.0,0.7).y());

    // srfOut.writeVRML("RR.wrl",Color(255,100,255),50,80);
}

TEST_F (object3DJoiningTest, testAsymRRAddAlongSeam){

    // Insert knots into the first surface
    Vector_FLOAT params(3);
    params[0] = 0.1;
    params[1] = 0.5;
    params[2] = 0.6;

    srf1A.refineKnotU(params);

    Object3D srfOut = map.joinSurfaces(srf1A,srf2A,"RR");

    // Test dimensions of output
    EXPECT_EQ(8,srfOut.ctrlPnts().rows());
    EXPECT_EQ(7,srfOut.ctrlPnts().cols());
    EXPECT_EQ(12,srfOut.knotU().size());
    EXPECT_EQ(11,srfOut.knotV().size());

}

TEST_F (object3DJoiningTest, testAsymDD){

    // Insert knots into the first surface
    Vector_FLOAT params(3);
    params[0] = 0.1;
    params[1] = 0.5;
    params[2] = 0.6;

    srf1D.refineKnotU(params);

    Object3D srfOut = map.joinSurfaces(srf1D,srf2A,"DD");

    // Test dimensions of output
    EXPECT_EQ(n_ctrl,srfOut.ctrlPnts().cols());
    EXPECT_EQ(2*n_ctrl - deg + 3,srfOut.ctrlPnts().rows());
    EXPECT_EQ(n_ctrl + 1 + deg,srfOut.knotV().size());
    EXPECT_EQ(2*n_ctrl + 1 + 3,srfOut.knotU().size());

    // Test Corner points 
    for (int i = 0; i < 4; i++){
        EXPECT_FLOAT_EQ(HPoint3Df(1.0,0.0,0.0,1.0).data[i],srfOut(0.0,0.0).data[i]);
        EXPECT_FLOAT_EQ(HPoint3Df(1.0,2.0,0.0,1.0).data[i],srfOut(1.0,0.0).data[i]);
        EXPECT_FLOAT_EQ(HPoint3Df(2.0,0.0,0.0,1.0).data[i],srfOut(0.0,1.0).data[i]);
        EXPECT_FLOAT_EQ(HPoint3Df(2.0,2.0,0.0,1.0).data[i],srfOut(1.0,1.0).data[i]);
    }
    

    // Test mid-extension points 
    EXPECT_FLOAT_EQ(1.0,srfOut(0.3,0.0).x());
    EXPECT_FLOAT_EQ(1.0,srfOut(0.7,0.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.3,1.0).x());
    EXPECT_FLOAT_EQ(2.0,srfOut(0.7,1.0).x());

    // srfOut.writeVRML("RR.wrl",Color(255,100,255),50,80);
}

TEST_F (object3DJoiningTest, testAsymDDAddAlongSeam){

    // Insert knots into the first surface
    Vector_FLOAT params(3);
    params[0] = 0.1;
    params[1] = 0.5;
    params[2] = 0.6;

    srf1D.refineKnotV(params);

    Object3D srfOut = map.joinSurfaces(srf1D,srf2A,"DD");

    // Test dimensions of output
    EXPECT_EQ(7,srfOut.ctrlPnts().rows());
    EXPECT_EQ(8,srfOut.ctrlPnts().cols());
    EXPECT_EQ(11,srfOut.knotU().size());
    EXPECT_EQ(12,srfOut.knotV().size());

}

//----------------------------------------------------------------
//----------------------------------------------------------------
//------------- MESH FROM SCAN TESTS  ------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
class meshFromScanTest : public ::testing::Test{
    protected:
        meshFromScanTest() : cloud(new pcl::PointCloud<pcl::PointNormal>(10,10)), nanArray(10,10), expectNanArray(10,10),
                rowFlags(10,1), colFlags(1,10), cloud_nan(new pcl::PointCloud<pcl::PointNormal>(10,10))
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
                    cloud_nan->at(j,i).x = NAN;
                    cloud_nan->at(j,i).y = NAN;
                    cloud_nan->at(j,i).z = NAN;
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

        pcl::PointCloud<pcl::PointNormal>::Ptr cloud;
        pcl::PointCloud<pcl::PointNormal>::Ptr cloud_nan;

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

    pcl::PointNormal expectP;

    expectP.x = 4.0f/9.0f;
    expectP.y = 5.0f/9.0f;
    expectP.z = 0.0f;

    cout << "Expected point:\n" << expectP << endl;
    
    cout << "Computed point:\n" << cloud->at(4,5) << endl;

    EXPECT_NEAR(expectP.x,cloud->at(4,5).x,0.0001);
    EXPECT_NEAR(expectP.y,cloud->at(4,5).y,0.0001);
    EXPECT_NEAR(expectP.z,cloud->at(4,5).z,0.0001);
    
    
}

TEST_F (meshFromScanTest, meshFromScanFunction){

    mp.numRowsDesired = 6;
    mp.numColsDesired = 6;

    pcl::PointCloud<pcl::PointNormal>::Ptr cloudOut (new pcl::PointCloud<pcl::PointNormal>(mp.numColsDesired,mp.numRowsDesired));

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

TEST_F (meshFromScanTest, tooFewRows){

    mp.numRowsDesired = 8;
    mp.numColsDesired = 8;

    mp.removeNanBuffer = 3;

    pcl::PointCloud<pcl::PointNormal>::Ptr cloudOut (new pcl::PointCloud<pcl::PointNormal>(mp.numColsDesired,mp.numRowsDesired));

    // Mesh from scan
    mp.meshFromScan(cloudOut, cloud);

    cout << "Cloud out size: height:" << cloudOut->height << ", width: " << cloudOut->width << endl;

    EXPECT_EQ(6,cloudOut->height);
    EXPECT_EQ(8,cloudOut->width);
    

}

TEST_F (meshFromScanTest, tooFewCols){

    mp.numRowsDesired = 6;
    mp.numColsDesired = 9;

    mp.removeNanBuffer = 3;

    pcl::PointCloud<pcl::PointNormal>::Ptr cloudOut (new pcl::PointCloud<pcl::PointNormal>(mp.numColsDesired,mp.numRowsDesired));

    // Mesh from scan
    mp.meshFromScan(cloudOut, cloud);

    cout << "Cloud out size: height:" << cloudOut->height << ", width: " << cloudOut->width << endl;

    EXPECT_EQ(6,cloudOut->height);
    EXPECT_EQ(7,cloudOut->width);

}

TEST_F (meshFromScanTest, allNans){

    pcl::PointCloud<pcl::PointNormal>::Ptr cloudOut (new pcl::PointCloud<pcl::PointNormal>(mp.numColsDesired,mp.numRowsDesired));

    // Mesh from scan
    mp.meshFromScan(cloudOut, cloud_nan);

    EXPECT_TRUE(mp.bRejectScan);

}

class meshFromScanLargeTest : public ::testing::Test{
    protected:
        meshFromScanLargeTest() : cloud(new pcl::PointCloud<pcl::PointNormal>(1000,1000))
        {
        
            int n_points = 1000;
            bool nanFlag = false;

            // Create data with Nans
            for (int i = 0; i<n_points; i++){
                for (int j = 0; j<n_points; j++){
                    if (i <= 1 || (i == 2 && (j < 3 || j > 7)) || i == 9 || j == 0 || i == 5 || j == 578 || i == 345 || j == 9 || (i == 5 && j == 4)){
                        nanFlag = true;
                    }
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

        pcl::PointCloud<pcl::PointNormal>::Ptr cloud;
  
        // Init class
        Mapping3D mp;
};

TEST_F (meshFromScanLargeTest, meshFromScanFunctionLarge){

    pcl::PointCloud<pcl::PointNormal>::Ptr cloudOut (new pcl::PointCloud<pcl::PointNormal>(mp.numColsDesired,mp.numRowsDesired));

    // Mesh from scan
    mp.meshFromScan(cloudOut, cloud);

    EXPECT_EQ(100,cloudOut->width);
    EXPECT_EQ(100,cloudOut->height);
}

//----------------------------------------------------------------
//----------------------------------------------------------------
//------------- Scan Processing Functions Test  ------------------
//----------------------------------------------------------------
//----------------------------------------------------------------

class scanProcessFunctionsTest : public ::testing::Test {
 protected:
  scanProcessFunctionsTest() : data(new pcl::PointCloud<pcl::PointNormal>(10,10))//, data_rot(new pcl::PointCloud<pcl::PointNormal>(10,10)) 
  {
    int n_points = 10;

    // *data = pcl::PointCloud<pcl::PointNormal>(10,10);
    
    // Create data
    for (int i = 0; i<n_points; i++){
        for (int j = 0; j<n_points; j++){
            data->at(j,i).x = (float)j/(n_points-1);
            data->at(j,i).y = (float)i/(n_points-1);
            data->at(j,i).z = 0.0f;
        } 
    }
  }

  pcl::PointCloud<pcl::PointNormal>::Ptr data;

  Mapping3D mp;
};

TEST_F (scanProcessFunctionsTest,computSearchMetricsTest){

    std::vector<float> searchMetrics = mp.computeSearchMetrics(data);

    EXPECT_FLOAT_EQ(0.5, searchMetrics[0]);
    EXPECT_FLOAT_EQ(0.5, searchMetrics[1]);
    EXPECT_FLOAT_EQ(0.0, searchMetrics[2]);
    EXPECT_FLOAT_EQ(1.0, searchMetrics[3]);
}

TEST_F (scanProcessFunctionsTest, addObjectTest){

    std::vector<float> searchMetrics = mp.computeSearchMetrics(data);

    mp.addObject(data,searchMetrics);

    EXPECT_EQ(1,mp.objectMap.size());
    EXPECT_EQ(1,mp.objectMetrics.size());

    EXPECT_FLOAT_EQ(0.0,mp.objectMap[0].pointAt(0.0,0.0).x());
    EXPECT_FLOAT_EQ(0.0,mp.objectMap[0].pointAt(0.0,0.0).y());
    EXPECT_FLOAT_EQ(0.0,mp.objectMap[0].pointAt(0.0,0.0).z());
    EXPECT_FLOAT_EQ(1.0,mp.objectMap[0].pointAt(1.0,1.0).x());
    EXPECT_FLOAT_EQ(1.0,mp.objectMap[0].pointAt(1.0,1.0).y());
    EXPECT_FLOAT_EQ(0.0,mp.objectMap[0].pointAt(1.0,1.0).z());


    // TDBM more tests?
}


//----------------------------------------------------------------
//----------------------------------------------------------------
//------------- Surface update testing ------------------
//----------------------------------------------------------------
//----------------------------------------------------------------

class updateSurfaceTest : public ::testing::Test {
 protected:
  updateSurfaceTest() : data(new pcl::PointCloud<pcl::PointNormal>(30,30)), data_rot(new pcl::PointCloud<pcl::PointNormal>(30,30)) 
  {
    
    int n_points = 30;

    // *data = pcl::PointCloud<pcl::PointNormal>(30,30);
    
    // Create data
    for (int i = 0; i<n_points; i++){
        for (int j = 0; j<n_points; j++){
            data->at(j,i).x = (float)j/(n_points-1);
            data->at(j,i).y = (float)i/(n_points-1);
            data->at(j,i).z = 0.0f;

            data_rot->at(j,i).x = (float)j/(n_points-1) - 0.5;
            data_rot->at(j,i).y = (float)i/(n_points-1) - 0.5;
            data_rot->at(j,i).z = 0.0f;
        } 
    }

    // Create transformation matrix
    Eigen::Affine3f transform = Eigen::Affine3f::Identity();

    // Compute Metrics
    std::vector<float> searchMetrics = mp.computeSearchMetrics(data);

    // Add object to map
    mp.addObject(data, searchMetrics);

    obj = mp.objectMap[0];

    cout << "data point after adding object is: " << data->at(0,0) << endl;

    // Transform data to get test data sets   
    // Left 
    transform.translation() << -0.4, 0.2, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_l[0], transform);
    transform.translation() << -0.4, -0.2, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_l[1], transform);
    transform.translation() << -0.4, 0.0, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_l[2], transform);  
    
    // Right
    transform.translation() << 0.4, 0.2, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_r[0], transform);
    transform.translation() << 0.4, -0.2, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_r[1], transform);
    transform.translation() << 0.4, 0.0, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_r[2], transform);

    // Down
    transform.translation() << 0.2, -0.4, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_d[0], transform);
    transform.translation() << -0.2, -0.4, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_d[1], transform);
    transform.translation() << 0.0, -0.4, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_d[2], transform);

    // Up
    transform.translation() << 0.2, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_u[0], transform);
    transform.translation() << -0.2, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_u[1], transform);
    transform.translation() << 0.0, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_u[2], transform);

    // Overlap (complete overlap)
    transform.scale(0.5f);// Scale by 1/2
    // cout << "transform for overlap is:\n" << transform.matrix() << endl;
    transform.translation() << 0.2, 0.2, 0.0;
    data_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_ol[0], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_ol[1], transform);
    transform.translation() << 0.4, 0.4, 0.0;
    data_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_ol[2], transform);

    // Mostly Overlap
    transform.translation() << -0.05, 0.2, 0.0;
    data_m_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_m_ol[0], transform);
    transform.translation() << 0.2, -0.05, 0.0;
    data_m_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_m_ol[1], transform);
    transform.translation() << 0.55, 0.55, 0.0;
    data_m_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_m_ol[2], transform);  

    // All new
    transform.scale(2.0f);// Scale back to 1.0
    transform.translation() << 1.5, 0.2, 0.0;
    data_all_new.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_all_new[0], transform);
    transform.translation() << 0.2, 1.5, 0.0;
    data_all_new.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_all_new[1], transform);
    transform.translation() << 1.5, 2.2, 0.0;
    data_all_new.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_all_new[2], transform);  

    // Rotations
    Eigen::Vector3f rotVec(0.0, 0.0, 1); // Rotations about the z axis
    // transform.scale(2.0f);// no scaling

    // Rotate 90 degrees (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (-M_PI/2.0f, rotVec));
    // cout << "transform for rot 90 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.4, 0.2, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot1[0], transform);
    transform.translation() << 0.4, 0.2, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot1[1], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot1[2], transform);
    transform.translation() << 0.2, -0.4, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot1[3], transform);
    
    // Rotate 90 degrees   (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (M_PI, rotVec));// Rotate 180 from 90
    // cout << "transform for rot -90 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.4, 0.2, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot2[0], transform);
    transform.translation() << 0.4, 0.2, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot2[1], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot2[2], transform);
    transform.translation() << 0.2, -0.4, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot2[3], transform);

    // Rotate 180 degrees   (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (M_PI/2, rotVec));// Rotate another 90 to get to 180
    // cout << "transform for rot 180 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.4, 0.2, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot3[0], transform);
    transform.translation() << 0.4, 0.2, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot3[1], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot3[2], transform);
    transform.translation() << 0.2, -0.4, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data_rot, *data_rot3[3], transform);
    
  }

  pcl::PointCloud<pcl::PointNormal>::Ptr data;
  pcl::PointCloud<pcl::PointNormal>::Ptr data_rot;



  // Vectors of point clouds 
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_l;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_r;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_d;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_u;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_ol;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_m_ol;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_all_new;
//   // For rotation
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot1;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot2;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot3;  

  // Class instatiation
  SplitSurface ss;
  Mapping3D mp;
  Object3D obj;
  Object3D obj_rot;

};

TEST_F (updateSurfaceTest, testLeft){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    for (int i = 0; i < 3; i++){
        cout << "Data point for data_l is: " << data_l[i]->at(0,0) << endl;

        mp.updateObject(0, data_l[i]);

        // Test dimensions of output
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
        
        // EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());

        // Test Corner points 
        // EXPECT_EQ(HPoint3Df(-0.4,0.2,0.0,1.0),mp.objectMap[0](0.0,0.0));
        // EXPECT_EQ(HPoint3Df(-0.4,1.2,0.0,1.0),mp.objectMap[0](1.0,0.0));
        // EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),mp.objectMap[0](0.0,1.0));
        // EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),mp.objectMap[0](1.0,1.0));
        switch (i){
            case 0 : 
                expectRes << -0.4, -0.4, 1.0, 1.0, 0.2, 1.2, 0.0, 1.0;
                break;
            case 1 : 
                expectRes << -0.4, -0.4, 1.0, 1.0, -0.2, 0.8, 0.0, 1.0;
                break;
            case 2 : 
                expectRes << -0.4, -0.4, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());

        // Reset object
        mp.updateObjectInMap(0, obj);

        // mp.objectMap[0].writeVRML("updateTest.wrl",Color(255,100,255),50,80);
    }
}

TEST_F (updateSurfaceTest, testRight){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    for (int i = 0; i < 3; i++){
        cout << "Data point for data_r is: " << data_r[i]->at(0,0) << endl;

        mp.updateObject(0, data_r[i]);

        // Test dimensions of output
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
        
        // EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());

        // Test Corner points 
        switch (i){
            case 0 : 
                expectRes << 0.0, 0.0, 1.4, 1.4, 0.0, 1.0, 0.2, 1.2;
                break;
            case 1 : 
                expectRes << 0.0, 0.0, 1.4, 1.4, 0.0, 1.0, -0.2, 0.8;
                break;
            case 2 : 
                expectRes << 0.0, 0.0, 1.4, 1.4, 0.0, 1.0, 0.0, 1.0;
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());

        // if (i == 1 ){
        //     mp.objectMap[0].writeVRML("updateTestR.wrl",Color(255,100,255),50,80);
        // }

        // Reset object
        mp.updateObjectInMap(0, obj);

        // mp.objectMap[0].writeVRML("updateTest.wrl",Color(255,100,255),50,80);
    }
}

TEST_F (updateSurfaceTest, testDown){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    for (int i = 0; i < 3; i++){
        cout << "Data point for data_d is: " << data_d[i]->at(0,0) << endl;

        mp.updateObject(0, data_d[i]);

        // Test dimensions of output
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
        
        // EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());

        // Test Corner points 
        // EXPECT_EQ(HPoint3Df(-0.4,0.2,0.0,1.0),mp.objectMap[0](0.0,0.0));
        // EXPECT_EQ(HPoint3Df(-0.4,1.2,0.0,1.0),mp.objectMap[0](1.0,0.0));
        // EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),mp.objectMap[0](0.0,1.0));
        // EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),mp.objectMap[0](1.0,1.0));
        switch (i){
            case 0 : 
                expectRes << 0.2, 0.0, 1.2, 1.0, -0.4, 1.0, -0.4, 1.0;
                break;
            case 1 : 
                expectRes << -0.2, 0.0, 0.8, 1.0, -0.4, 1.0, -0.4, 1.0;
                break;
            case 2 : 
                expectRes << 0.0, 0.0, 1.0, 1.0, -0.4, 1.0, -0.4, 1.0;
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());

        // if (i == 1 ){
        //     mp.objectMap[0].writeVRML("updateTestd.wrl",Color(255,100,255),50,80);
        // }

        // Reset object
        mp.updateObjectInMap(0, obj);

    }
}

TEST_F (updateSurfaceTest, testUp){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    for (int i = 0; i < 3; i++){
        cout << "Data point for data_u is: " << data_u[i]->at(0,0) << endl;

        mp.updateObject(0, data_u[i]);

        // Test dimensions of output
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
        
        // EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());

        // Test Corner points 
        // EXPECT_EQ(HPoint3Df(-0.4,0.2,0.0,1.0),mp.objectMap[0](0.0,0.0));
        // EXPECT_EQ(HPoint3Df(-0.4,1.2,0.0,1.0),mp.objectMap[0](1.0,0.0));
        // EXPECT_EQ(HPoint3Df(1.0,0.0,0.0,1.0),mp.objectMap[0](0.0,1.0));
        // EXPECT_EQ(HPoint3Df(1.0,1.0,0.0,1.0),mp.objectMap[0](1.0,1.0));
        switch (i){
            case 0 : 
                expectRes << 0.0, 0.2, 1.0, 1.2, 0.0, 1.4, 0.0, 1.4;
                break;
            case 1 : 
                expectRes << 0.0, -0.2, 1.0, 0.8, 0.0, 1.4, 0.0, 1.4;
                break;
            case 2 : 
                expectRes << 0.0, 0.0, 1.0, 1.0, 0.0, 1.4, 0.0, 1.4;
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());

        if (i == 1 ){
            mp.objectMap[0].writeVRML("updateTestU.wrl",Color(255,100,255),50,80);
        }
        

        // Reset object
        mp.updateObjectInMap(0, obj);
    }
}

TEST_F (updateSurfaceTest, testOverlap){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    for (int i = 0; i < 3; i++){
        cout << "Data point for data_ol is: " << data_ol[i]->at(0,0) << endl;

        mp.updateObject(0, data_ol[i]);

        // Test dimensions of output
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
        

        expectRes << 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0;

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());


        // Reset object
        mp.updateObjectInMap(0, obj);
    }
}

TEST_F (updateSurfaceTest, testMostlyOverlap){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    for (int i = 0; i < 3; i++){
        cout << "Data point for data_m_ol is: " << data_m_ol[i]->at(0,0) << endl;

        mp.updateObject(0, data_m_ol[i]);

        // Test dimensions of output
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
        EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
        

        expectRes << 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0;

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());


        // Reset object
        mp.updateObjectInMap(0, obj);
    }
}

TEST_F (updateSurfaceTest, testAllNew){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    for (int i = 0; i < 3; i++){
        cout << "Data point for data_all_new is: " << data_all_new[i]->at(0,0) << endl;

        mp.updateObject(0, data_all_new[i]);

        
        
        switch (i){
            case 0 : 
                expectRes << 0.0, 0.0, 2.5, 2.5, 0.0, 1.0, 0.2, 1.2;
                // Test dimensions of output
                EXPECT_EQ(30,mp.objectMap[0].ctrlPnts().cols());
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
                break;
            case 1 : 
                expectRes << 0.0, 0.2, 1.0, 1.2, 0.0, 2.5, 0.0, 2.5;
                // Test dimensions of output
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                EXPECT_EQ(30,mp.objectMap[0].ctrlPnts().rows());
                break;
            case 2 : 
                expectRes << 0.0, 1.5, 1.0, 2.5, 0.0, 3.2, 0.0, 3.2;
                // Test dimensions of output
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                EXPECT_EQ(30,mp.objectMap[0].ctrlPnts().rows());
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());


        // Reset object
        mp.updateObjectInMap(0, obj);
    }
}

TEST_F (updateSurfaceTest, testRotatePlus90){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    // Compute Metrics
    std::vector<float> searchMetrics = mp.computeSearchMetrics(data_rot);

    mp.objectMap.erase(mp.objectMap.begin());
    mp.objectMetrics.erase(mp.objectMetrics.begin());

    // Add object to map
    mp.addObject(data_rot, searchMetrics);

    obj_rot = mp.objectMap[0];

    for (int i = 0; i < 4; i++){
        cout << "Data point for data_rot1 is: " << data_rot1[i]->at(0,0) << endl;

        

        mp.updateObject(0, data_rot1[i]);

        // Test dimensions of output
        
        
        // EXPECT_EQ(2*n_ctrl - deg,srfOut.ctrlPnts().cols());


        switch (i){
            case 0 : 
                // LD
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
                expectRes << -0.9, -0.9, 0.5, 0.5, -0.3, 0.7, -0.5, 0.5;
                break;
            case 1 : 
                // RU
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
                expectRes << -0.5, -0.5, 0.9, 0.9, -0.5, 0.5, -0.3, 0.7;
                break;
            case 2 : 
                // UL
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                expectRes << -0.5, -0.3, 0.5, 0.7, -0.5, 0.9, -0.5, 0.9;
                break;
            case 3 : 
                // DR
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                expectRes << -0.3, -0.5, 0.7, 0.5, -0.9, 0.5, -0.9, 0.5;
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());
       
        mp.updateObjectInMap(0, obj_rot);
        
    }
}

TEST_F (updateSurfaceTest, testRotateMinus90){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    // Compute Metrics
    std::vector<float> searchMetrics = mp.computeSearchMetrics(data_rot);

    mp.objectMap.erase(mp.objectMap.begin());
    mp.objectMetrics.erase(mp.objectMetrics.begin());

    // Add object to map
    mp.addObject(data_rot, searchMetrics);

    obj_rot = mp.objectMap[0];

    for (int i = 0; i < 4; i++){
        cout << "Data point for data_rot2 is: " << data_rot2[i]->at(0,0) << endl;

        mp.updateObject(0, data_rot2[i]);

        // Test dimensions of output

        switch (i){
            case 0 : 
                // LD
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
                expectRes << -0.9, -0.9, 0.5, 0.5, -0.3, 0.7, -0.5, 0.5;
                break;
            case 1 : 
                // RU
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
                expectRes << -0.5, -0.5, 0.9, 0.9, -0.5, 0.5, -0.3, 0.7;
                break;
            case 2 : 
                // UL
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                expectRes << -0.5, -0.3, 0.5, 0.7, -0.5, 0.9, -0.5, 0.9;
                break;
            case 3 : 
                // DR
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                expectRes << -0.3, -0.5, 0.7, 0.5, -0.9, 0.5, -0.9, 0.5;
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());
       

        // Reset object
        mp.updateObjectInMap(0, obj_rot);
    }
}

TEST_F (updateSurfaceTest, testRotate180){

    Eigen::Array<float, 1, 8> expectRes(1,8);

    // Compute Metrics
    std::vector<float> searchMetrics = mp.computeSearchMetrics(data_rot);

    mp.objectMap.erase(mp.objectMap.begin());
    mp.objectMetrics.erase(mp.objectMetrics.begin());

    // Add object to map
    mp.addObject(data_rot, searchMetrics);

    obj_rot = mp.objectMap[0];

    for (int i = 0; i < 4; i++){
        cout << "Data point for data_rot3 is: " << data_rot3[i]->at(0,0) << endl;

        mp.updateObject(0, data_rot3[i]);

        // Test dimensions of output

        switch (i){
            case 0 : 
                // LD
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
                expectRes << -0.9, -0.9, 0.5, 0.5, -0.3, 0.7, -0.5, 0.5;
                break;
            case 1 : 
                // RU
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().rows());
                expectRes << -0.5, -0.5, 0.9, 0.9, -0.5, 0.5, -0.3, 0.7;
                break;
            case 2 : 
                // UL
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                expectRes << -0.5, -0.3, 0.5, 0.7, -0.5, 0.9, -0.5, 0.9;
                break;
            case 3 : 
                // DR
                EXPECT_EQ(15,mp.objectMap[0].ctrlPnts().cols());
                expectRes << -0.3, -0.5, 0.7, 0.5, -0.9, 0.5, -0.9, 0.5;
                break;
        }

        EXPECT_FLOAT_EQ(expectRes[0],mp.objectMap[0](0.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[1],mp.objectMap[0](1.0,0.0).x());
        EXPECT_FLOAT_EQ(expectRes[2],mp.objectMap[0](0.0,1.0).x());
        EXPECT_FLOAT_EQ(expectRes[3],mp.objectMap[0](1.0,1.0).x());

        EXPECT_FLOAT_EQ(expectRes[4],mp.objectMap[0](0.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[5],mp.objectMap[0](1.0,0.0).y());
        EXPECT_FLOAT_EQ(expectRes[6],mp.objectMap[0](0.0,1.0).y());
        EXPECT_FLOAT_EQ(expectRes[7],mp.objectMap[0](1.0,1.0).y());
       

        // Reset object
        mp.updateObjectInMap(0, obj_rot);
    }
}


//----------------------------------------------------------------
//----------------------------------------------------------------
//------------- Surface update testing ------------------
//----------------------------------------------------------------
//----------------------------------------------------------------

class nurbsDataFromPCNonRect : public ::testing::Test {
 protected:
  nurbsDataFromPCNonRect() : data(new pcl::PointCloud<pcl::PointNormal>(30,30)),
   data_rot(new pcl::PointCloud<pcl::PointNormal>(30,30)),
   dataC(new pcl::PointCloud<pcl::PointNormal>(30,30)),
   dataC2(new pcl::PointCloud<pcl::PointNormal>(30,30))
  {
    
    int n_points = 30;
    float offset;
    float offset2;

    // *data = pcl::PointCloud<pcl::PointNormal>(30,30);
    
    // Create data
    for (int i = 0; i<n_points; i++){
        if (i < 15){
          offset = i*0.3/14.0;
        }else{
          offset = (n_points - i - 1)*0.3/14.0;
        }
        // cout << "Offset is: " << offset << endl;
        for (int j = 0; j<n_points; j++){
            data->at(j,i).x = (float)j/(n_points-1);
            data->at(j,i).y = (float)i/(n_points-1);
            data->at(j,i).z = 0.0f;

            dataC->at(j,i).x = (float)j/(n_points-1) - offset;
            dataC->at(j,i).y = (float)i/(n_points-1);
            dataC->at(j,i).z = 0.0f;

            if (j < 15){
                offset2 = j*0.3/14.0;
            }else{
                offset2 = (n_points - j - 1)*0.3/14.0;
            }

            dataC2->at(j,i).x = (float)j/(n_points-1);
            dataC2->at(j,i).y = (float)i/(n_points-1) - offset2;
            dataC2->at(j,i).z = 0.0f;

            data_rot->at(j,i).x = (float)j/(n_points-1) - 0.5;
            data_rot->at(j,i).y = (float)i/(n_points-1) - 0.5;
            data_rot->at(j,i).z = 0.0f;
        } 
    }

    // Create transformation matrix
    Eigen::Affine3f transform = Eigen::Affine3f::Identity();

    // Compute Metrics
    std::vector<float> searchMetrics = mp.computeSearchMetrics(data);

    // Add object to map
    mp.addObject(dataC, searchMetrics);

    // Compute Metrics
    searchMetrics = mp.computeSearchMetrics(dataC);

    // Add object to map
    mp.addObject(dataC, searchMetrics);
    searchMetrics = mp.computeSearchMetrics(dataC2);
    mp.addObject(dataC2, searchMetrics);

    obj = mp.objectMap[0];
    objC = mp.objectMap[1];
    objC2 = mp.objectMap[2];
    

    cout << "data point after adding object is: " << data->at(0,0) << endl;

    // Transform data to get test data sets   
    // Left 
    transform.translation() << -0.7, 0.2, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_l[0], transform);
    transform.translation() << -0.7, -0.2, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_l[1], transform);
    transform.translation() << -0.7, 0.0, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_l[2], transform);
    transform.translation() << -1.1, 0.0, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_l[3], transform);
    
    // Right
    transform.translation() << 0.4, 0.2, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_r[0], transform);
    transform.translation() << 0.4, -0.2, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_r[1], transform);
    transform.translation() << 0.4, 0.0, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_r[2], transform);
    transform.translation() << 0.8, 0.0, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_r[3], transform);
    
    


    // Down
    transform.translation() << 0.2, -0.7, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_d[0], transform);
    transform.translation() << -0.2, -0.7, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_d[1], transform);
    transform.translation() << 0.0, -0.7, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_d[2], transform);

    // Up
    transform.translation() << 0.2, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_u[0], transform);
    transform.translation() << -0.2, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_u[1], transform);
    transform.translation() << 0.0, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(30,30)));
    pcl::transformPointCloud(*data, *data_u[2], transform);

    // Rotations
    Eigen::Vector3f rotVec(0.0, 0.0, 1); // Rotations about the z axis
    transform.scale(1.0f);// no scaling

    // Rotate 90 degrees (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (-M_PI/2.0f, rotVec));
    // cout << "transform for rot 90 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.7, 0.2, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[0], transform);
    transform.translation() << 0.9, 0.2, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[1], transform);
    transform.translation() << 0.2, 0.9, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[2], transform);
    transform.translation() << 0.2, -0.7, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[3], transform);
    
    // Rotate 90 degrees   (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (M_PI, rotVec));// Rotate 180 from 90
    // cout << "transform for rot -90 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.7, 0.2, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[0], transform);
    transform.translation() << 0.9, 0.2, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[1], transform);
    transform.translation() << 0.2, 0.9, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[2], transform);
    transform.translation() << 0.2, -0.7, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[3], transform);

    // Rotate 180 degrees   (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (M_PI/2, rotVec));// Rotate another 90 to get to 180
    // cout << "transform for rot 180 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.7, 0.2, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot3[0], transform);
    transform.translation() << 0.9, 0.2, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot3[1], transform);
    transform.translation() << 0.2, 0.9, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot3[2], transform);
    transform.translation() << 0.2, -0.7, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot3[3], transform);


  }

  pcl::PointCloud<pcl::PointNormal>::Ptr data;
  pcl::PointCloud<pcl::PointNormal>::Ptr dataC;
  pcl::PointCloud<pcl::PointNormal>::Ptr dataC2;
  pcl::PointCloud<pcl::PointNormal>::Ptr data_rot;



  // Vectors of point clouds 
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_l;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_r;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_d;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_u;
  
  // For rotation
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot1;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot2;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot3;  

  // Class instatiation
  SplitSurface ss;
  Mapping3D mp;
  Object3D obj;
  Object3D objC;
  Object3D objC2;

};

TEST_F (nurbsDataFromPCNonRect, testNormalRight){

    int ms, mt;
    int row1,row2;

    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface right
        ss.splitNewSurfaceObservation(data,data_r[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_r[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        switch (j){
            case 0: row1 = 0; row2 = 23; break;
            case 1: row1 = 6; row2 = 29; break;
            case 2: row1 = 0; row2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).x,mesh(0,mt-1).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).y,mesh(0,mt-1).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).z,mesh(0,mt-1).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).x,mesh(ms-1,mt-1).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).y,mesh(ms-1,mt-1).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).z,mesh(ms-1,mt-1).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(17,row1).x,mesh(0,0).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row1).y,mesh(0,0).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row1).z,mesh(0,0).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(17,row2).x,mesh(ms-1,0).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row2).y,mesh(ms-1,0).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row2).z,mesh(ms-1,0).z());
    }
}

TEST_F (nurbsDataFromPCNonRect, testNormalLeft){

    int ms, mt;
    int row1,row2;

    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface right
        ss.splitNewSurfaceObservation(data,data_l[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_l[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        switch (j){
            case 0: row1 = 0; row2 = 23; break;
            case 1: row1 = 6; row2 = 29; break;
            case 2: row1 = 0; row2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_l[j]->at(20,row1).x,mesh(0,mt-1).x());
        EXPECT_FLOAT_EQ(data_l[j]->at(20,row1).y,mesh(0,mt-1).y());
        EXPECT_FLOAT_EQ(data_l[j]->at(20,row1).z,mesh(0,mt-1).z());

        EXPECT_FLOAT_EQ(data_l[j]->at(20,row2).x,mesh(ms-1,mt-1).x());
        EXPECT_FLOAT_EQ(data_l[j]->at(20,row2).y,mesh(ms-1,mt-1).y());
        EXPECT_FLOAT_EQ(data_l[j]->at(20,row2).z,mesh(ms-1,mt-1).z());

        EXPECT_FLOAT_EQ(data_l[j]->at(0,row1).x,mesh(0,0).x());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row1).y,mesh(0,0).y());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row1).z,mesh(0,0).z());

        EXPECT_FLOAT_EQ(data_l[j]->at(0,row2).x,mesh(ms-1,0).x());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row2).y,mesh(ms-1,0).y());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row2).z,mesh(ms-1,0).z());
    }
}

TEST_F (nurbsDataFromPCNonRect, testNormalDown){

    int ms, mt;
    int col1,col2;

    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface 
        ss.splitNewSurfaceObservation(data,data_d[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_d[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        switch (j){
            case 0: col1 = 0; col2 = 23; break;
            case 1: col1 = 6; col2 = 29; break;
            case 2: col1 = 0; col2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_d[j]->at(col2,0).x,mesh(0,mt-1).x());
        EXPECT_FLOAT_EQ(data_d[j]->at(col2,0).y,mesh(0,mt-1).y());
        EXPECT_FLOAT_EQ(data_d[j]->at(col2,0).z,mesh(0,mt-1).z());
        
        EXPECT_FLOAT_EQ(data_d[j]->at(col2,20).x,mesh(ms-1,mt-1).x());
        EXPECT_FLOAT_EQ(data_d[j]->at(col2,20).y,mesh(ms-1,mt-1).y());
        EXPECT_FLOAT_EQ(data_d[j]->at(col2,20).z,mesh(ms-1,mt-1).z());

        EXPECT_FLOAT_EQ(data_d[j]->at(col1,0).x,mesh(0,0).x());
        EXPECT_FLOAT_EQ(data_d[j]->at(col1,0).y,mesh(0,0).y());
        EXPECT_FLOAT_EQ(data_d[j]->at(col1,0).z,mesh(0,0).z());

        EXPECT_FLOAT_EQ(data_d[j]->at(col1,20).x,mesh(ms-1,0).x());
        EXPECT_FLOAT_EQ(data_d[j]->at(col1,20).y,mesh(ms-1,0).y());
        EXPECT_FLOAT_EQ(data_d[j]->at(col1,20).z,mesh(ms-1,0).z());

    }
}

TEST_F (nurbsDataFromPCNonRect, testNormalUp){

    int ms, mt;
    int col1,col2;

    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface 
        ss.splitNewSurfaceObservation(data,data_u[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_u[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        switch (j){
            case 0: col1 = 0; col2 = 23; break;
            case 1: col1 = 6; col2 = 29; break;
            case 2: col1 = 0; col2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_u[j]->at(col2,17).x,mesh(0,mt-1).x());
        EXPECT_FLOAT_EQ(data_u[j]->at(col2,17).y,mesh(0,mt-1).y());
        EXPECT_FLOAT_EQ(data_u[j]->at(col2,17).z,mesh(0,mt-1).z());
        
        EXPECT_FLOAT_EQ(data_u[j]->at(col2,data->width-1).x,mesh(ms-1,mt-1).x());
        EXPECT_FLOAT_EQ(data_u[j]->at(col2,data->width-1).y,mesh(ms-1,mt-1).y());
        EXPECT_FLOAT_EQ(data_u[j]->at(col2,data->width-1).z,mesh(ms-1,mt-1).z());

        EXPECT_FLOAT_EQ(data_u[j]->at(col1,17).x,mesh(0,0).x());
        EXPECT_FLOAT_EQ(data_u[j]->at(col1,17).y,mesh(0,0).y());
        EXPECT_FLOAT_EQ(data_u[j]->at(col1,17).z,mesh(0,0).z());

        EXPECT_FLOAT_EQ(data_u[j]->at(col1,data->width-1).x,mesh(ms-1,0).x());
        EXPECT_FLOAT_EQ(data_u[j]->at(col1,data->width-1).y,mesh(ms-1,0).y());
        EXPECT_FLOAT_EQ(data_u[j]->at(col1,data->width-1).z,mesh(ms-1,0).z());

    }
}

TEST_F (nurbsDataFromPCNonRect, testRight){

    int ms, mt;
    int row1, row2, rowSelect;

    Eigen::Array<float,1,Eigen::Dynamic> expectInd(1,data->height);
    expectInd << 17,17,16,16,15,14,14,13,12,12,11,11,10,9,9,9,9,10,11,11,12,12,13,14,14,15,16,16,17,17;

    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface right
        ss.splitNewSurfaceObservation(dataC,data_r[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_r[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        // Select row1 and row 2 to give expected values
        switch (j){
            case 0: row1 = 0; row2 = 23; break;
            case 1: row1 = 6; row2 = 29; break;
            case 2: row1 = 0; row2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).x,mesh(0,mt-1).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).y,mesh(0,mt-1).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).z,mesh(0,mt-1).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).x,mesh(ms-1,mt-1).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).y,mesh(ms-1,mt-1).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).z,mesh(ms-1,mt-1).z());
        
        // int rowCount = 0;
        
        // for (int i = 0; i < ms; i++){
        //     if (j == 0 && i<6){
        //         rowSelect = row1;
        //     }else if (j == 1 && i > 22){
        //         rowSelect = row2;
        //     }else{
        //         rowSelect = row1 + rowCount;
        //         rowCount++;
        //     }
        //     cout << "Row Select is: " << rowSelect << endl;
        //     EXPECT_FLOAT_EQ(data_r[j]->at(expectInd[i],rowSelect).x,mesh(i,0).x());
        //     EXPECT_FLOAT_EQ(data_r[j]->at(expectInd[i],rowSelect).y,mesh(i,0).y());
        //     EXPECT_FLOAT_EQ(data_r[j]->at(expectInd[i],rowSelect).z,mesh(i,0).z());
        // }
        
    }
}

TEST_F (nurbsDataFromPCNonRect, testLeft){

    int ms, mt;
    int row1, row2, rowSelect;

    Eigen::Array<float,1,Eigen::Dynamic> expectInd(1,data->height);
    expectInd << 20,20,19,18,18,17,17,16,15,15,14,13,13,12,12,12,12,13,13,14,15,15,16,17,17,18,18,19,20,20;

    // Get the indices
    for (int j = 2; j < 3; j++){
        // Split surface right
        ss.splitNewSurfaceObservation(dataC,data_l[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_l[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        // Select row1 and row 2 to give expected values
        switch (j){
            case 0: row1 = 0; row2 = 23; break;
            case 1: row1 = 6; row2 = 29; break;
            case 2: row1 = 0; row2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_l[j]->at(0,row1).x,mesh(0,0).x());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row1).y,mesh(0,0).y());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row1).z,mesh(0,0).z());

        EXPECT_FLOAT_EQ(data_l[j]->at(0,row2).x,mesh(ms-1,0).x());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row2).y,mesh(ms-1,0).y());
        EXPECT_FLOAT_EQ(data_l[j]->at(0,row2).z,mesh(ms-1,0).z());
        
        // int rowCount = 0;
        
        // for (int i = 0; i < ms; i++){
        //     if (j == 0 && i<6){
        //         rowSelect = row1;
        //     }else if (j == 1 && i > 22){
        //         rowSelect = row2;
        //     }else{
        //         rowSelect = row1 + rowCount;
        //         rowCount++;
        //     }
        //     cout << "Row Select is: " << rowSelect << endl;
        //     EXPECT_FLOAT_EQ(data_l[j]->at(expectInd[i],rowSelect).x,mesh(i,mt-1).x());
        //     EXPECT_FLOAT_EQ(data_l[j]->at(expectInd[i],rowSelect).y,mesh(i,mt-1).y());
        //     EXPECT_FLOAT_EQ(data_l[j]->at(expectInd[i],rowSelect).z,mesh(i,mt-1).z());
        // }
        
    }
}

TEST_F (nurbsDataFromPCNonRect, testDown){

    int ms, mt;
    int col1,col2;
    int colSelect;

    Eigen::Array<float,1,Eigen::Dynamic> expectInd(1,data->width);
    expectInd << 20,20,19,18,18,17,17,16,15,15,14,13,13,12,12,12,12,13,13,14,15,15,16,17,17,18,18,19,20,20;


    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface 
        ss.splitNewSurfaceObservation(dataC2,data_d[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_d[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        switch (j){
            case 0: col1 = 0; col2 = 23; break;
            case 1: col1 = 6; col2 = 29; break;
            case 2: col1 = 0; col2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_d[j]->at(col2,0).x,mesh(0,mt-1).x());
        EXPECT_FLOAT_EQ(data_d[j]->at(col2,0).y,mesh(0,mt-1).y());
        EXPECT_FLOAT_EQ(data_d[j]->at(col2,0).z,mesh(0,mt-1).z());

        EXPECT_FLOAT_EQ(data_d[j]->at(col1,0).x,mesh(0,0).x());
        EXPECT_FLOAT_EQ(data_d[j]->at(col1,0).y,mesh(0,0).y());
        EXPECT_FLOAT_EQ(data_d[j]->at(col1,0).z,mesh(0,0).z());

        // int colCount = 0;       
        
        // for (int i = 0; i < mt; i++){
        //     if (j == 0 && i<6){
        //         colSelect = col1;
        //     }else if (j == 1 && i > 22){
        //         colSelect = col2;
        //     }else{
        //         colSelect = col1 + colCount;
        //         colCount++;
        //     }
        //     cout << "col Select is: " << colSelect << endl;
        //     EXPECT_FLOAT_EQ(data_d[j]->at(colSelect, expectInd[i]).x,mesh(ms-1,i).x());
        //     EXPECT_FLOAT_EQ(data_d[j]->at(colSelect, expectInd[i]).y,mesh(ms-1,i).y());
        //     EXPECT_FLOAT_EQ(data_d[j]->at(colSelect, expectInd[i]).z,mesh(ms-1,i).z());
        // }
    }
}

TEST_F (nurbsDataFromPCNonRect, testUp){

    int ms, mt;
    int col1,col2;
    int colSelect;

    Eigen::Array<float,1,Eigen::Dynamic> expectInd(1,data->width);
    expectInd << 17, 17, 16, 16, 15, 14, 14, 13, 12, 12, 11, 11, 10, 9, 9, 9, 9, 10, 11, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 17;


    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface 
        ss.splitNewSurfaceObservation(dataC2,data_u[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_u[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        switch (j){
            case 0: col1 = 0; col2 = 23; break;
            case 1: col1 = 6; col2 = 29; break;
            case 2: col1 = 0; col2 = data->width-1; break;
        }

        EXPECT_FLOAT_EQ(data_u[j]->at(col2,data->height-1).x,mesh(ms-1,mt-1).x());
        EXPECT_FLOAT_EQ(data_u[j]->at(col2,data->height-1).y,mesh(ms-1,mt-1).y());
        EXPECT_FLOAT_EQ(data_u[j]->at(col2,data->height-1).z,mesh(ms-1,mt-1).z());

        EXPECT_FLOAT_EQ(data_u[j]->at(col1,data->height-1).x,mesh(ms-1,0).x());
        EXPECT_FLOAT_EQ(data_u[j]->at(col1,data->height-1).y,mesh(ms-1,0).y());
        EXPECT_FLOAT_EQ(data_u[j]->at(col1,data->height-1).z,mesh(ms-1,0).z());

        // int colCount = 0;       
        
        // for (int i = 0; i < mt; i++){
        //     if (j == 0 && i<6){
        //         colSelect = col1;
        //     }else if (j == 1 && i > 22){
        //         colSelect = col2;
        //     }else{
        //         colSelect = col1 + colCount;
        //         colCount++;
        //     }
        //     cout << "col Select is: " << colSelect << endl;
        //     EXPECT_FLOAT_EQ(data_u[j]->at(colSelect, expectInd[i]).x,mesh(0,i).x());
        //     EXPECT_FLOAT_EQ(data_u[j]->at(colSelect, expectInd[i]).y,mesh(0,i).y());
        //     EXPECT_FLOAT_EQ(data_u[j]->at(colSelect, expectInd[i]).z,mesh(0,i).z());
        // }
    }
}

TEST_F (nurbsDataFromPCNonRect, testRightObj){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?
    int ms, mt;
    int row1, row2, rowSelect;

    // Get the indices
    for (int j = 0; j < 3; j++){
        // Split surface right
        ss.splitNewSurfaceObservation(dataC,data_r[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_r[j], ss.newRowIndices, ss.newColIndices);

        // Create object
        Object3D objNew(mesh);

        // // Write to file
        // if (j == 1){
        //     objNew.writeVRML("RnonRect1.wrl",Color(255,100,255),50,80);
        // }
    }
}

TEST_F (nurbsDataFromPCNonRect, testRightSurfUpdate){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?
    int row1,row2;
    mp.useNonRectData = true;

    pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>(45,45,pcl::PointNormal()));

    std::string filename = "/home/bjm/SpaceCRAFT/Results/nurbs_overlap_initial.pcd";
    mp.writeObjectPCDFile(filename.c_str(), 1, 10, 10);

    pcl::PCDWriter writer;
    filename = "/home/bjm/SpaceCRAFT/Results/nurbs_overlap_original_data.pcd";
    writer.write<pcl::PointNormal> (filename, *dataC, false);
    filename = "/home/bjm/SpaceCRAFT/Results/nurbs_overlap_new_data.pcd";
    writer.write<pcl::PointNormal> (filename, *data_r[0], false);

    // Get the indices
    for (int j = 0; j < 1; j++){
        
        // UPdate object 2
        mp.updateObject(1, data_r[j]);


        switch (j){
            case 0: row1 = 0; row2 = 23; break;
            case 1: row1 = 6; row2 = 29; break;
            case 2: row1 = 0; row2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).x,mp.objectMap[1](0.0,1.0).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).y,mp.objectMap[1](0.0,1.0).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).z,mp.objectMap[1](0.0,1.0).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).x,mp.objectMap[1](1.0,1.0).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).y,mp.objectMap[1](1.0,1.0).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).z,mp.objectMap[1](1.0,1.0).z());

        EXPECT_FLOAT_EQ(data->at(0,0).x,mp.objectMap[1](0.0,0.0).x());
        EXPECT_FLOAT_EQ(data->at(0,0).y,mp.objectMap[1](0.0,0.0).y());
        EXPECT_FLOAT_EQ(data->at(0,0).z,mp.objectMap[1](0.0,0.0).z());

        EXPECT_FLOAT_EQ(data->at(0,data->height-1).x,mp.objectMap[1](1.0,0.0).x());
        EXPECT_FLOAT_EQ(data->at(0,data->height-1).y,mp.objectMap[1](1.0,0.0).y());
        EXPECT_FLOAT_EQ(data->at(0,data->height-1).z,mp.objectMap[1](1.0,0.0).z());
        
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_r[j]);
        Object3D objd(mesh);


        EXPECT_FLOAT_EQ(objd.ctrlPnts()(0,9).x(),mp.objectMap[1].ctrlPnts()(0,mp.objectMap[1].ctrlPnts().cols()-1).x());

        // Test generating data - to test it doesn't break. 
        mp.pointCloudFromObject3D(1,45,45,cloud);
        bool bIsnoNan = true;
        for (int i = 0; i < 45*45; i++){
            if (!pcl::isFinite(cloud->at(i))){
                bIsnoNan = false;
            }
        }
        EXPECT_TRUE(bIsnoNan);
        // Write to file
        if (j == 0){
            mp.objectMap[1].writeVRML("UpdatedSurf0.wrl",Color(255,100,255),50,80);
            filename = "/home/bjm/SpaceCRAFT/Results/nurbs_overlap_final.pcd";
            mp.writeObjectPCDFile(filename.c_str(), 1, 10, 15);
        }else if (j == 1){
            mp.objectMap[1].writeVRML("UpdatedSurf1.wrl",Color(255,100,255),50,80);
        }else if (j == 2){
            filename = "/home/bjm/SpaceCRAFT/Results/nurbs_overlap_final.pcd";
            mp.writeObjectPCDFile(filename.c_str(), 1, 25, 55);
        }
        // Reset object
        mp.updateObjectInMap(1, objC);
    }
}
/*
TEST_F (nurbsDataFromPCNonRect, testLeftSurfUpdate){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 3; j++){
        
        // UPdate object 2
        mp.updateObject(1, data_l[j]);

        // Write to file
        // if (j == 0){
        //     mp.objectMap[1].writeVRML("UpdatedSurf0L.wrl",Color(255,100,255),50,80);
        // }else if (j == 1){
        //     mp.objectMap[1].writeVRML("UpdatedSurf1L.wrl",Color(255,100,255),50,80);
        // }else{
        //     mp.objectMap[1].writeVRML("UpdatedSurf2L.wrl",Color(255,100,255),50,80);
        // }
        // Reset object
        mp.updateObjectInMap(1, objC);
    }
}

TEST_F (nurbsDataFromPCNonRect, testDownSurfUpdate){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 3; j++){
        
        // UPdate object 2
        mp.updateObject(2, data_d[j]);

        // Write to file
        if (j == 0){
            mp.objectMap[2].writeVRML("UpdatedSurf0D.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[2].writeVRML("UpdatedSurf1D.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[2].writeVRML("UpdatedSurf2D.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(2, objC2);
    }
}

TEST_F (nurbsDataFromPCNonRect, testUpSurfUpdate){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 3; j++){
        
        // UPdate object 2
        mp.updateObject(2, data_u[j]);

        // Write to file
        if (j == 0){
            mp.objectMap[2].writeVRML("UpdatedSurf0U.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[2].writeVRML("UpdatedSurf1U.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[2].writeVRML("UpdatedSurf2U.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(2, objC2);
    }
}

TEST_F (nurbsDataFromPCNonRect, testSurfUpUpdateRotate90){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 4; j++){
        
        // UPdate object 2
        mp.updateObject(2, data_rot1[j]);

        // Write to file
        if (j == 0){
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot0.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot1.wrl",Color(255,100,255),50,80);
        }else if (j == 2){
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot2.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot3.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(2, objC2);
        
    }
}

TEST_F (nurbsDataFromPCNonRect, testSurfRightUpdateRotate90){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 4; j++){
        
        // UPdate object 2
        mp.updateObject(1, data_rot1[j]);

        // Write to file
        if (j == 0){
            mp.objectMap[1].writeVRML("UpdatedSurfR_rot0.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[1].writeVRML("UpdatedSurfR_rot1.wrl",Color(255,100,255),50,80);
        }else if (j == 2){
            mp.objectMap[1].writeVRML("UpdatedSurfR_rot2.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[1].writeVRML("UpdatedSurfR_rot3.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(1, objC);
        
    }
}

TEST_F (nurbsDataFromPCNonRect, testUpSurfUpdateRotateMinus90){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 4; j++){
        
        // UPdate object 2
        mp.updateObject(2, data_rot1[j]);

        // Write to file
        if (j == 0){
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot0.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot1.wrl",Color(255,100,255),50,80);
        }else if (j == 2){
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot2.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[2].writeVRML("UpdatedSurfU_rot3.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(2, objC2);
        
    }
}

TEST_F (nurbsDataFromPCNonRect, testUpSurfUpdateDown){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 1; j++){
        
        // UPdate object 2
        mp.updateObject(2, data_rot3[2]);

        // Write to file
        if (j == 0){
            mp.objectMap[2].writeVRML("UpdatedSurf0UD.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[2].writeVRML("UpdatedSurf1UD.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[2].writeVRML("UpdatedSurf2UD.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(2, objC2);
        
    }
}

TEST_F (nurbsDataFromPCNonRect, testUpSurfUpdateRight){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 1; j++){
        
        // UPdate object 2
        mp.updateObject(2, data_rot2[2]);

        // Write to file
        if (j == 0){
            mp.objectMap[2].writeVRML("UpdatedSurf0UR.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[2].writeVRML("UpdatedSurf1UR.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[2].writeVRML("UpdatedSurf2UR.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(2, objC2);
    }
}

TEST_F (nurbsDataFromPCNonRect, testUpSurfUpdateLeft){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;

    // Get the indices
    for (int j = 0; j < 1; j++){
        
        // UPdate object 2
        mp.updateObject(2, data_rot1[2]);

        // Write to file
        if (j == 0){
            mp.objectMap[2].writeVRML("UpdatedSurf0UL.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[2].writeVRML("UpdatedSurf1UL.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[2].writeVRML("UpdatedSurf2UL.wrl",Color(255,100,255),50,80);
        }
        // Reset object
        mp.updateObjectInMap(2, objC2);
    }
}

TEST_F (nurbsDataFromPCNonRect, testRightSurfUpdateSuccessiveCompact){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;
    int nPoints = 45;

    pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>(nPoints,nPoints,pcl::PointNormal()));


    // Get the indices
    for (int j = 0; j < 3; j++){
        
        // UPdate object 2
        mp.updateObject(1, data_r[j]);

        // Update again
        mp.updateObject(1, data_r[3]);
        
        // Test generating data - to test it doesn't break. 
        mp.pointCloudFromObject3D(1,nPoints,nPoints,cloud);
        bool bIsnoNan = true;
        for (int i = 0; i < nPoints*nPoints; i++){
            if (!pcl::isFinite(cloud->at(i))){
                bIsnoNan = false;
                // cout << "\n\n\t\tNANS IN SUCCESSIVE SURFACE!!\n\n" << endl;
                cout << i << endl;
            }
        }
        EXPECT_TRUE(bIsnoNan);

        // Write to file
        if (j == 0){
            mp.objectMap[1].writeVRML("UpdatedSurf0.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[1].writeVRML("UpdatedSurf1.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[1].writeVRML("UpdatedSurf2.wrl",Color(255,100,255),50,80);
        }

        // Reset object
        mp.updateObjectInMap(1, objC);
    }
}

TEST_F (nurbsDataFromPCNonRect, testRightSurfUpdateSuccessiveDetails){
    // Just to make sure it doesn't crash
    // TODO -some tests on the output object?

    mp.useNonRectData = true;
    int nPoints = 45;

    pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>(nPoints,nPoints,pcl::PointNormal()));


    // Get the indices
    for (int j = 0; j < 1; j++){
        
        // UPdate object 2
        mp.updateObject(1, data_r[j]);

        // Update again
        
        cout << "Getting point cloud" << endl;
        pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC(new pcl::PointCloud<pcl::PointNormal>(125, 125, pcl::PointNormal()));
        mp.pointCloudFromObject3D(1, 125, 125, mapObjPC);

        cout << "Splitting surface " << endl;
        // Split surface right
        ss.splitNewSurfaceObservation(mapObjPC,data_r[3]);

        cout << "Getting data indices" << endl;
        // Get new data indices 
        ss.getNewDataIndices();
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        cout << "getting mesh from PC" << endl;

        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_r[3], ss.newRowIndices, ss.newColIndices);
        
        cout << "Checking Mesh, with rows, cols: " << mesh.rows() << ", " << mesh.cols() << endl;
        cout << mesh << endl;
        for (int i = 0; i < mesh.rows(); i++){
            for (int j = 0; j < mesh.cols(); j++){
                // cout << "in mesh check lop ij is " << i << ", " << j << endl;
                if (!std::isfinite(mesh(i,j).x()) || !std::isfinite(mesh(i,j).y()) || !std::isfinite(mesh(i,j).z())){
                    cout << "Mesh nan at (i,j): (" << i << ", " << j << ").\n";
                }
            }
        }
        cout << "Done checking mesh " << endl;

        EXPECT_FLOAT_EQ(data_r[3]->at(29,6).x,mesh(0,mesh.cols()-1).x());
        EXPECT_FLOAT_EQ(data_r[3]->at(29,6).y,mesh(0,mesh.cols()-1).y());
        EXPECT_FLOAT_EQ(data_r[3]->at(29,6).z,mesh(0,mesh.cols()-1).z());

        EXPECT_FLOAT_EQ(data_r[3]->at(29,29).x,mesh(mesh.rows()-1,mesh.cols()-1).x());
        EXPECT_FLOAT_EQ(data_r[3]->at(29,29).y,mesh(mesh.rows()-1,mesh.cols()-1).y());
        EXPECT_FLOAT_EQ(data_r[3]->at(29,29).z,mesh(mesh.rows()-1,mesh.cols()-1).z());

        // Get point cloud
        pcl::PointCloud<pcl::PointNormal>::Ptr newCloud(new pcl::PointCloud<pcl::PointNormal>(mesh.cols(), mesh.rows(), pcl::PointNormal()));

        for (int i = 0; i < mesh.rows(); i++){
            for (int j = 0 ; j < mesh.cols(); j++){
                newCloud->at(j,i).x = mesh(i,j).x();
                newCloud->at(j,i).y = mesh(i,j).y();
                newCloud->at(j,i).z = mesh(i,j).z();
            }
        }
        // Write cloud
        pcl::PCDWriter writer;
        writer.write<pcl::PointNormal> ("testNewPartError.pcd", *newCloud, false);

        

        std::vector<int> nCtrlNew = mp.computeNumberOfControlPoints(ss.extendDirection, mesh, mp.objectMap[1]);

        Object3D obj(mesh, 3, 3, nCtrlNew[0], nCtrlNew[1]);

        EXPECT_FLOAT_EQ(obj(0.0,1.0).x(),mesh(0,mesh.cols()-1).x());
        EXPECT_FLOAT_EQ(obj(0.0,1.0).y(),mesh(0,mesh.cols()-1).y());
        EXPECT_FLOAT_EQ(obj(0.0,1.0).z(),mesh(0,mesh.cols()-1).z());

        EXPECT_FLOAT_EQ(obj(0.0,0.0).x(),mesh(0,0).x());
        EXPECT_FLOAT_EQ(obj(0.0,0.0).y(),mesh(0,0).y());

        EXPECT_FLOAT_EQ(obj(1.0,1.0).x(),mesh(mesh.rows()-1,mesh.cols()-1).x());
        EXPECT_FLOAT_EQ(obj(1.0,1.0).y(),mesh(mesh.rows()-1,mesh.cols()-1).y());
        EXPECT_FLOAT_EQ(obj(1.0,1.0).z(),mesh(mesh.rows()-1,mesh.cols()-1).z());

        EXPECT_FLOAT_EQ(obj(1.0,0.0).x(),mesh(mesh.rows()-1,0).x());
        EXPECT_FLOAT_EQ(obj(1.0,0.0).y(),mesh(mesh.rows()-1,0).y());


        // Test generating data - to test it doesn't break. 
        mp.pointCloudFromObject3D(1,nPoints,nPoints,cloud);
        bool bIsnoNan = true;
        for (int i = 0; i < nPoints*nPoints; i++){
            if (!pcl::isFinite(cloud->at(i))){
                bIsnoNan = false;
                // cout << "\n\n\t\tNANS IN SUCCESSIVE SURFACE!!\n\n" << endl;
                cout << i << endl;
            }
        }
        EXPECT_TRUE(bIsnoNan);

        // Write to file
        if (j == 0){
            mp.objectMap[1].writeVRML("UpdatedSurf0.wrl",Color(255,100,255),50,80);
        }else if (j == 1){
            mp.objectMap[1].writeVRML("UpdatedSurf1.wrl",Color(255,100,255),50,80);
        }else{
            mp.objectMap[1].writeVRML("UpdatedSurf2.wrl",Color(255,100,255),50,80);
        }

        // Reset object
        mp.updateObjectInMap(1, objC);
    }
}


TEST_F (nurbsDataFromPCNonRect, testMultipleIndices){

    int ms, mt;
    int row1,row2;

    // Get the indices
    for (int j = 2; j < 3; j++){
        // Split surface right
        ss.splitNewSurfaceObservation(data,data_r[j]);

        // Get new data indices 
        ss.getNewDataIndices();

        Eigen::Array<int,1,Eigen::Dynamic> indices(1,ss.newRowIndices.rows());
        for (int i = 0; i < ss.newRowIndices.rows(); i++){
            if (i < 13){
                ss.newRowIndices(i,0) = 0;
                ss.newRowIndices(i,1) = 0;
            }else if (i < 20){
                ss.newRowIndices(i,0) = 1;
                ss.newRowIndices(i,1) = 1;
            }
        }

        // Print
        cout << "New Row indices are:\n" << ss.newRowIndices << endl;
        cout << "New Col indices are:\n" << ss.newColIndices << endl;

        // Get mesh out 
        Matrix_Point3Df mesh = mp.nurbsDataFromPointCloud(data_r[j], ss.newRowIndices, ss.newColIndices);

        ms = mesh.rows();
        mt = mesh.cols();

        cout << "In test ms, mt are: " << ms << ", " << mt << endl;

        switch (j){
            case 0: row1 = 0; row2 = 23; break;
            case 1: row1 = 6; row2 = 29; break;
            case 2: row1 = 0; row2 = data->height-1; break;
        }

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).x,mesh(0,mt-1).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).y,mesh(0,mt-1).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row1).z,mesh(0,mt-1).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).x,mesh(ms-1,mt-1).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).y,mesh(ms-1,mt-1).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(data->width-1,row2).z,mesh(ms-1,mt-1).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(17,row1).x,mesh(0,0).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row1).y,mesh(0,0).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row1).z,mesh(0,0).z());

        EXPECT_FLOAT_EQ(data_r[j]->at(17,row2).x,mesh(ms-1,0).x());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row2).y,mesh(ms-1,0).y());
        EXPECT_FLOAT_EQ(data_r[j]->at(17,row2).z,mesh(ms-1,0).z());

    }
}

*/
} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
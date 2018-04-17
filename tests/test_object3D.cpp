
#include "cans/object3D.h"
#include "gtest/gtest.h"


namespace {

using namespace PLib;

class object3DTest : public ::testing::Test {
 protected:
  object3DTest() : scan(3,3), scan2(25,25) {
    float k;
  
    for (int i = 0; i<scan.rows(); i++){
      for (int j=0; j<scan.cols(); j++){
        k = i*j + 2.0*i - j*3.0;
        scan(i,j) = Point3Df(i,j,k);
      } 
    }

    k = 0;
    
    for (int i = 0; i<scan2.rows(); i++){
        for (int j=0; j<scan2.cols(); j++){
        scan2(i,j) = Point3Df(i,j,k);
        }
    }

    // Least Squares Fit
    surf.leastSquares(scan2,3,3,10,10) ;
  }

  ~object3DTest() {;}


  Matrix_Point3Df scan;
  Matrix_Point3Df scan2;
  PlNurbsSurfacef surf ;

};

// Constructors

TEST_F (object3DTest, BlankConstructor){
    Object3D obj;
}

TEST_F (object3DTest, DataConstructorDefault){
    Object3D obj(scan2);

    // Test a few points to make sure the surface fitting is the same
    EXPECT_EQ(obj.pointAt(0.1,0.1),surf.pointAt(0.1,0.1));
    EXPECT_EQ(obj.pointAt(0.1,0.5),surf.pointAt(0.1,0.5));
    EXPECT_EQ(obj.pointAt(0.3,0.7),surf.pointAt(0.3,0.7));
}


TEST_F (object3DTest, DataConstructor){
    Object3D obj(scan2,3,3,7,7);

    // Least Squares Fit
    surf.leastSquares(scan2,3,3,7,7) ;

    // Test a few points to make sure the surface fitting is the same
    EXPECT_EQ(obj.pointAt(0.1,0.1),surf.pointAt(0.1,0.1));
    EXPECT_EQ(obj.pointAt(0.1,0.5),surf.pointAt(0.1,0.5));
    EXPECT_EQ(obj.pointAt(0.3,0.7),surf.pointAt(0.3,0.7));
}

TEST_F (object3DTest, readAndWrite){

    // Create the object
    Object3D obj(scan2);

    // Write it to file
    obj.write("test_filename.obj");

    // Create a second object
    Object3D obj2;

    // Read it 
    obj2.readObject3D("test_filename.obj");
    
    // Compare centres
    EXPECT_NEAR(obj.getObjSize(),obj2.getObjSize(),0.1);
    EXPECT_NEAR(obj.getCentre().x(),obj2.getCentre().x(),1.0);
    EXPECT_NEAR(obj.getCentre().y(),obj2.getCentre().y(),1.0);
    EXPECT_NEAR(obj.getCentre().z(),obj2.getCentre().z(),0.1);
}

// Getters
TEST_F (object3DTest, getCentre){
    Object3D obj;

    Point3Df expect(0.0,0.0,0.0);
    Point3Df centre = obj.getCentre();

    EXPECT_EQ(centre,expect);
}

TEST_F (object3DTest, getColor){
    Object3D obj;

    Point3Df expect(0.0,0.0,0.0);
    Point3Df color = obj.getColor();

    EXPECT_EQ(color,expect);
}

TEST_F (object3DTest, getObjSize){
    Object3D obj;

    float expect(0.0);
    float obj_size = obj.getObjSize();

    EXPECT_EQ(obj_size,expect);
}


TEST_F (object3DTest, TestComputeCentreDataSimple){
    // Initialise class
    Object3D obj(scan);

    // Get centre
    Point3Df centre = obj.getCentre();
    Point3Df expect(1,1,0);

    EXPECT_EQ(centre, expect);
}

TEST_F (object3DTest, TestComputeCentreData){
    // Initialise class
    Object3D obj(scan2);

    // Get centre
    Point3Df centre = obj.getCentre();
    Point3Df expect(12.0,12.0,0);

    EXPECT_EQ(centre, expect);
}

TEST_F (object3DTest, TestComputeCentreCP){
    // Initialise class
    Object3D obj(scan2);

    // compute centre with the control points
    obj.computeCentreFromControlPoints();

    // Get centre
    Point3Df centre = obj.getCentre();
    // Point3Df expect(11.1,11.1,0);

    // EXPECT_EQ(centre, expect);
}

TEST_F (object3DTest, computeSizeData){
    // Initialise class
    Object3D obj(scan2);

    // Get size 
    float expect(24.0);
    float obj_size = obj.getObjSize();

    EXPECT_EQ(obj_size,expect);
}

TEST_F (object3DTest, computeSizeDataCP){
    // Initialise class
    Object3D obj(scan2);

    // Get size 
    float expect(24.0);
    obj.computeSizeFromControlPoints();
    float obj_size = obj.getObjSize();

    // EXPECT_FLOAT_EQ(obj_size,expect); // Note that EXPECT_FLOAT_EQ won't work for user defined types
    EXPECT_NEAR(obj_size,expect,1e-4); // Note that EXPECT_FLOAT_EQ won't work for user defined types
}


} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
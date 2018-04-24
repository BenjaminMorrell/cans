
#include "cans/object3D.h"
#include "gtest/gtest.h"


namespace {

using namespace PLib;

class object3DTest : public ::testing::Test {
 protected:
  object3DTest() : scan(3,3), scan2(25,25), scan3(25,25) {
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

    k = 2;

    for (int i = 0; i<scan2.rows(); i++){
        for (int j=0; j<scan3.cols(); j++){
            scan3(i,j) = Point3Df(i,j,k);
        }
    }


    // Least Squares Fit
    surf.leastSquares(scan2,3,3,10,10) ;
  }

  ~object3DTest() {;}


  Matrix_Point3Df scan;
  Matrix_Point3Df scan2;
  Matrix_Point3Df scan3;
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

TEST_F (object3DTest, KnotsConstructor){
    Object3D obj(scan2,3,3,7,7);

    Vector_FLOAT knotU = obj.knotU();
    Vector_FLOAT knotV = obj.knotV();
    Matrix_HPoint3Df controlP = obj.ctrlPnts();

    // Use knot and control points constructor
    Object3D obj2(obj.degreeU(),obj.degreeV(),knotU,knotV,controlP);

    // Test a few points to make sure the surface fitting is the same
    EXPECT_EQ(obj2.pointAt(0.1,0.1),obj.pointAt(0.1,0.1));
    EXPECT_EQ(obj2.pointAt(0.1,0.5),obj.pointAt(0.1,0.5));
    EXPECT_EQ(obj2.pointAt(0.3,0.7),obj.pointAt(0.3,0.7));
}

TEST_F (object3DTest, KnotsConstructorEigen){
    Object3D obj(scan2,3,3,7,7);

    int sizeKnots = obj.knotU().size();

    // Form Eigen structures
    Eigen::Array<float,1,Eigen::Dynamic> Uvec(1, sizeKnots);
    Eigen::Array<float,1,Eigen::Dynamic> Vvec(1, sizeKnots);
    Eigen::Array<float,3,Eigen::Dynamic> controlP(3,obj.ctrlPnts().rows()*obj.ctrlPnts().cols());

    // Fill Eigen arrays with values. 
    for (int i = 0; i < sizeKnots; i++){
        Uvec[i] = obj.knotU()[i];
        Vvec[i] = obj.knotV()[i];
    }

    int index = 0;
    for (int i = 0; i < obj.ctrlPnts().rows(); i++){
        for (int j = 0; j < obj.ctrlPnts().cols(); j++){
            controlP(0,index) = obj.ctrlPnts()(i,j).x();
            controlP(1,index) = obj.ctrlPnts()(i,j).y();
            controlP(2,index) = obj.ctrlPnts()(i,j).z();
            index++;
        }
    }

    // cout << "Control points are: " << controlP << endl;

    // cout << "Control points shoulld be:\n" << obj.ctrlPnts() << endl;


    // Use knot and control points constructor
    Object3D obj2(obj.degreeU(),obj.degreeV(), Uvec, Vvec, controlP, obj.ctrlPnts().rows(),obj.ctrlPnts().cols());

    // Test a few points to make sure the surface fitting is the same
    EXPECT_EQ(obj2.pointAt(0.1,0.1),obj.pointAt(0.1,0.1));
    EXPECT_EQ(obj2.pointAt(0.1,0.5),obj.pointAt(0.1,0.5));
    EXPECT_EQ(obj2.pointAt(0.3,0.7),obj.pointAt(0.3,0.7));
}

TEST_F (object3DTest, KnotsConstructorEigenAsym){
    Object3D obj(scan2,3,3,7,10);

    int sizeKnotU = obj.knotU().size();
    int sizeKnotV = obj.knotV().size();

    // Form Eigen structures
    Eigen::Array<float,1,Eigen::Dynamic> Uvec(1, sizeKnotU);
    Eigen::Array<float,1,Eigen::Dynamic> Vvec(1, sizeKnotV);
    Eigen::Array<float,3,Eigen::Dynamic> controlP(3,obj.ctrlPnts().rows()*obj.ctrlPnts().cols());

    // Fill Eigen arrays with values. 
    for (int i = 0; i < sizeKnotU; i++){
        Uvec[i] = obj.knotU()[i];
    }
    for (int i = 0; i < sizeKnotV; i++){
        Vvec[i] = obj.knotV()[i];
    }

    int index = 0;
    for (int i = 0; i < obj.ctrlPnts().rows(); i++){
        for (int j = 0; j < obj.ctrlPnts().cols(); j++){
            controlP(0,index) = obj.ctrlPnts()(i,j).x();
            controlP(1,index) = obj.ctrlPnts()(i,j).y();
            controlP(2,index) = obj.ctrlPnts()(i,j).z();
            index++;
        }
    }

    // cout << "Control points are: " << controlP << endl;

    // cout << "Control points shoulld be:\n" << obj.ctrlPnts() << endl;


    // Use knot and control points constructor
    Object3D obj2(obj.degreeU(),obj.degreeV(), Uvec, Vvec, controlP, obj.ctrlPnts().rows(),obj.ctrlPnts().cols());

    // Test a few points to make sure the surface fitting is the same
    EXPECT_EQ(obj2.pointAt(0.1,0.1),obj.pointAt(0.1,0.1));
    EXPECT_EQ(obj2.pointAt(0.1,0.5),obj.pointAt(0.1,0.5));
    EXPECT_EQ(obj2.pointAt(0.3,0.7),obj.pointAt(0.3,0.7));
}

TEST_F (object3DTest, updateObject3D){

    // Create 2 objects
    Object3D obj1(scan2);
    Object3D obj2(scan3);

    // Create knot vector to use
    int nKnotU = obj1.knotU().size();
    Eigen::Array<float,1,Eigen::Dynamic> knotU(1,nKnotU);
    int nKnotV = obj1.knotV().size();
    Eigen::Array<float,1,Eigen::Dynamic> knotV(1,nKnotV);

    int nCtrlS = obj1.ctrlPnts().rows();
    int nCtrlT = obj1.ctrlPnts().cols();

    // Create control vector to use
    Eigen::Array<float,3,Eigen::Dynamic> controlP(3,nCtrlS*nCtrlT);

    // Fill information from obj1
    for (int i = 0; i < nKnotU; i++){
        knotU(i) = obj1.knotU()[i];
    }
    for (int i = 0; i < nKnotV; i++){
        knotV(i) = obj1.knotV()[i];
    }
    int index = 0;
    for (int i = 0; i < nCtrlS; i++){
        for (int j = 0; j < nCtrlT; j++){
            controlP(0,index) = obj1.ctrlPnts()(i,j).x();
            controlP(1,index) = obj1.ctrlPnts()(i,j).y();
            controlP(2,index) = obj1.ctrlPnts()(i,j).z();
            index++;
        }
    }
    // cout << "End index is " << index << ", and number of control points is: " << nCtrlS*nCtrlT << endl;

    // Update obj2 with parameters of obj 1
    obj2.updateObject3D(obj1.degreeU(),obj1.degreeV(),knotU,knotV,controlP,nCtrlS,nCtrlT);

    // Test result
    // Test some points
    EXPECT_EQ(obj2.pointAt(0.1,0.1),obj1.pointAt(0.1,0.1));
    EXPECT_EQ(obj2.pointAt(0.1,0.5),obj1.pointAt(0.1,0.5));
    EXPECT_EQ(obj2.pointAt(0.3,0.7),obj1.pointAt(0.3,0.7));

    // Check control points
    EXPECT_EQ(obj2.ctrlPnts()(3,5),obj1.ctrlPnts()(3,5));
    EXPECT_EQ(obj2.ctrlPnts()(0,0),obj1.ctrlPnts()(0,0));
    EXPECT_EQ(obj2.ctrlPnts()(6,7),obj1.ctrlPnts()(6,7));

    // Check knots
    EXPECT_EQ(obj2.knotU()[5],obj1.knotU()[5]);
    EXPECT_EQ(obj2.knotU()[7],obj1.knotU()[7]);
    EXPECT_EQ(obj2.knotV()[5],obj1.knotV()[5]);
    EXPECT_EQ(obj2.knotV()[7],obj1.knotV()[7]);
}

TEST_F (object3DTest, updateObject3DStdVec){

    // Create 2 objects
    Object3D obj1(scan2);
    Object3D obj2(scan3);

    std::vector<float> knotU;
    std::vector<float> knotV;
    for (int i = 0; i < obj1.knotU().size(); i++){
        knotU.push_back(obj1.knotU()[i]);
    }
    for (int i = 0; i < obj1.knotV().size(); i++){
        knotV.push_back(obj1.knotV()[i]);
    }
    int nCtrlS = obj1.ctrlPnts().rows();
    int nCtrlT = obj1.ctrlPnts().cols();

    std::vector<float> controlX;
    std::vector<float> controlY;
    std::vector<float> controlZ;

    for (int i = 0; i < nCtrlS; i++){
        for (int j = 0; j < nCtrlT; j++){
            controlX.push_back(obj1.ctrlPnts()(i,j).x());
            controlY.push_back(obj1.ctrlPnts()(i,j).y());
            controlZ.push_back(obj1.ctrlPnts()(i,j).z());
        }
    }
    // cout << "End index is " << index << ", and number of control points is: " << nCtrlS*nCtrlT << endl;

    // Update obj2 with parameters of obj 1
    obj2.updateObject3DCPP(obj1.degreeU(),obj1.degreeV(),knotU,knotV,controlX,controlY,controlZ,nCtrlS,nCtrlT);

    // Test result
    // Test some points
    EXPECT_EQ(obj2.pointAt(0.1,0.1),obj1.pointAt(0.1,0.1));
    EXPECT_EQ(obj2.pointAt(0.1,0.5),obj1.pointAt(0.1,0.5));
    EXPECT_EQ(obj2.pointAt(0.3,0.7),obj1.pointAt(0.3,0.7));

    // Check control points
    EXPECT_EQ(obj2.ctrlPnts()(3,5),obj1.ctrlPnts()(3,5));
    EXPECT_EQ(obj2.ctrlPnts()(0,0),obj1.ctrlPnts()(0,0));
    EXPECT_EQ(obj2.ctrlPnts()(6,7),obj1.ctrlPnts()(6,7));

    // Check knots
    EXPECT_EQ(obj2.knotU()[5],obj1.knotU()[5]);
    EXPECT_EQ(obj2.knotU()[7],obj1.knotU()[7]);
    EXPECT_EQ(obj2.knotV()[5],obj1.knotV()[5]);
    EXPECT_EQ(obj2.knotV()[7],obj1.knotV()[7]);
}

TEST_F (object3DTest, updateObject3DAsymObj){

    // Create 2 objects
    Object3D obj1(scan2, 3, 3, 7, 10);
    Object3D obj2(scan3);

    // Create knot vector to use
    int nKnotU = obj1.knotU().size();
    Eigen::Array<float,1,Eigen::Dynamic> knotU(1,nKnotU);
    int nKnotV = obj1.knotV().size();
    Eigen::Array<float,1,Eigen::Dynamic> knotV(1,nKnotV);

    int nCtrlS = obj1.ctrlPnts().rows();
    int nCtrlT = obj1.ctrlPnts().cols();

    // Create control vector to use
    Eigen::Array<float,3,Eigen::Dynamic> controlP(3,nCtrlS*nCtrlT);

    // Fill information from obj1
    for (int i = 0; i < nKnotU; i++){
        knotU(i) = obj1.knotU()[i];
    }
    for (int i = 0; i < nKnotV; i++){
        knotV(i) = obj1.knotV()[i];
    }
    int index = 0;
    for (int i = 0; i < nCtrlS; i++){
        for (int j = 0; j < nCtrlT; j++){
            controlP(0,index) = obj1.ctrlPnts()(i,j).x();
            controlP(1,index) = obj1.ctrlPnts()(i,j).y();
            controlP(2,index) = obj1.ctrlPnts()(i,j).z();
            index++;
        }
    }
    cout << "End index is " << index << ", and number of control points is: " << nCtrlS*nCtrlT << endl;

    // Update obj2 with parameters of obj 1
    obj2.updateObject3D(obj1.degreeU(),obj1.degreeV(),knotU,knotV,controlP,nCtrlS,nCtrlT);

    cout << "updated the object " << endl;
    // Test result
    // Test some points
    EXPECT_EQ(obj2.pointAt(0.1,0.1),obj1.pointAt(0.1,0.1));
    EXPECT_EQ(obj2.pointAt(0.1,0.5),obj1.pointAt(0.1,0.5));
    EXPECT_EQ(obj2.pointAt(0.3,0.7),obj1.pointAt(0.3,0.7));

    // Check control points
    EXPECT_EQ(obj2.ctrlPnts()(3,5),obj1.ctrlPnts()(3,5));
    EXPECT_EQ(obj2.ctrlPnts()(0,0),obj1.ctrlPnts()(0,0));
    EXPECT_EQ(obj2.ctrlPnts()(6,7),obj1.ctrlPnts()(6,7));

    // Check knots
    EXPECT_EQ(obj2.knotU()[5],obj1.knotU()[5]);
    EXPECT_EQ(obj2.knotU()[7],obj1.knotU()[7]);
    EXPECT_EQ(obj2.knotV()[5],obj1.knotV()[5]);
    EXPECT_EQ(obj2.knotV()[7],obj1.knotV()[7]);
}

TEST_F (object3DTest, updateObject3DStdVecAsym){

    // Create 2 objects
    Object3D obj1(scan2, 3, 3, 7, 10);
    Object3D obj2(scan3);

    std::vector<float> knotU;
    std::vector<float> knotV;
    for (int i = 0; i < obj1.knotU().size(); i++){
        knotU.push_back(obj1.knotU()[i]);
    }
    for (int i = 0; i < obj1.knotV().size(); i++){
        knotV.push_back(obj1.knotV()[i]);
    }
    int nCtrlS = obj1.ctrlPnts().rows();
    int nCtrlT = obj1.ctrlPnts().cols();

    std::vector<float> controlX;
    std::vector<float> controlY;
    std::vector<float> controlZ;

    for (int i = 0; i < nCtrlS; i++){
        for (int j = 0; j < nCtrlT; j++){
            controlX.push_back(obj1.ctrlPnts()(i,j).x());
            controlY.push_back(obj1.ctrlPnts()(i,j).y());
            controlZ.push_back(obj1.ctrlPnts()(i,j).z());
        }
    }
    // cout << "End index is " << index << ", and number of control points is: " << nCtrlS*nCtrlT << endl;

    // Update obj2 with parameters of obj 1
    obj2.updateObject3DCPP(obj1.degreeU(),obj1.degreeV(),knotU,knotV,controlX,controlY,controlZ,nCtrlS,nCtrlT);

    // Test result
    // Test some points
    EXPECT_EQ(obj2.pointAt(0.1,0.1),obj1.pointAt(0.1,0.1));
    EXPECT_EQ(obj2.pointAt(0.1,0.5),obj1.pointAt(0.1,0.5));
    EXPECT_EQ(obj2.pointAt(0.3,0.7),obj1.pointAt(0.3,0.7));

    // Check control points
    EXPECT_EQ(obj2.ctrlPnts()(3,5),obj1.ctrlPnts()(3,5));
    EXPECT_EQ(obj2.ctrlPnts()(0,0),obj1.ctrlPnts()(0,0));
    EXPECT_EQ(obj2.ctrlPnts()(6,7),obj1.ctrlPnts()(6,7));

    // Check knots
    EXPECT_EQ(obj2.knotU()[5],obj1.knotU()[5]);
    EXPECT_EQ(obj2.knotU()[7],obj1.knotU()[7]);
    EXPECT_EQ(obj2.knotV()[5],obj1.knotV()[5]);
    EXPECT_EQ(obj2.knotV()[7],obj1.knotV()[7]);
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

TEST_F (object3DTest, getSurfacePointCloud){
    // Initialise class
    Object3D obj(scan2);
    
    // Number of points
    int msSurf = 25;
    int mtSurf = 25;

    // Init a point cloud
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>(mtSurf,msSurf,pcl::PointNormal())); 

    // Get the point cloud
    obj.getSurfacePointCloud(cloud,msSurf,mtSurf);

    float val(0.0001);

    // Test corner points
    EXPECT_FLOAT_EQ(cloud->at(0,0).x,obj.pointAt(val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,0).y,obj.pointAt(val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,0).z,obj.pointAt(val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).x,obj.pointAt(1.0-val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).y,obj.pointAt(1.0-val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).z,obj.pointAt(1.0-val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).x,obj.pointAt(val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).y,obj.pointAt(val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).z,obj.pointAt(val,1.0-val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).x,obj.pointAt(1.0-val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).y,obj.pointAt(1.0-val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).z,obj.pointAt(1.0-val,1.0-val).z());

}

TEST_F (object3DTest, getSurfacePointCloudSkewedMsMt){
    // Initialise class
    Object3D obj(scan2);
    
    // Number of points
    int msSurf = 45;
    int mtSurf = 25;

    // Init a point cloud
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>(mtSurf,msSurf,pcl::PointNormal())); 

    // Get the point cloud
    obj.getSurfacePointCloud(cloud,msSurf,mtSurf);

    float val(0.0001);

    // Test corner points
    EXPECT_FLOAT_EQ(cloud->at(0,0).x,obj.pointAt(val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,0).y,obj.pointAt(val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,0).z,obj.pointAt(val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).x,obj.pointAt(1.0-val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).y,obj.pointAt(1.0-val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).z,obj.pointAt(1.0-val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).x,obj.pointAt(val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).y,obj.pointAt(val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).z,obj.pointAt(val,1.0-val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).x,obj.pointAt(1.0-val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).y,obj.pointAt(1.0-val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).z,obj.pointAt(1.0-val,1.0-val).z());

}

TEST_F (object3DTest, getSurfacePointCloudAsymmetricalObject){
    // Initialise class
    Object3D obj(scan2, 3, 3, 7, 10);
    
    // Number of points
    int msSurf = 25;
    int mtSurf = 25;

    // Init a point cloud
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>(mtSurf,msSurf,pcl::PointNormal())); 

    // Get the point cloud
    obj.getSurfacePointCloud(cloud,msSurf,mtSurf);

    float val(0.0001);

    // Test corner points
    EXPECT_FLOAT_EQ(cloud->at(0,0).x,obj.pointAt(val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,0).y,obj.pointAt(val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,0).z,obj.pointAt(val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).x,obj.pointAt(1.0-val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).y,obj.pointAt(1.0-val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).z,obj.pointAt(1.0-val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).x,obj.pointAt(val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).y,obj.pointAt(val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).z,obj.pointAt(val,1.0-val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).x,obj.pointAt(1.0-val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).y,obj.pointAt(1.0-val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).z,obj.pointAt(1.0-val,1.0-val).z());

}


TEST_F (object3DTest, getDistanceToPointDefault){
    // Initialise class
    Object3D obj(scan2);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(2.0,-dist,0.05);

}

TEST_F (object3DTest, getDistanceToPoint){
    // Initialise class
    Object3D obj(scan2);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query,55,55);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(-2.0,dist,0.1);

}

TEST_F (object3DTest, getDistanceToPointSkewedMsMt){
    // Initialise class
    Object3D obj(scan2);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query,35,55);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(-2.0,dist,0.1);

}

TEST_F (object3DTest, getDistanceToPointAsymmetricObj){
    // Initialise class
    Object3D obj(scan2, 3, 3, 7, 10);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query,35,35);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(-2.0,dist,0.1);

}

TEST_F (object3DTest, getDistanceToPointAsym2){
    // Initialise class
    Object3D obj(scan2,3,3,5,12);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    cout << "\n\nObj size is: (" << obj.ctrlPnts().rows() << ", " << obj.ctrlPnts().cols() << ")\n\n";


    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(-2.0,dist,0.1);

}

TEST_F (object3DTest, getBatchDistanceToSurface){
    // Initialise class
    Object3D obj(scan2);

    // Test batch of points
    int nPoints = 35;
    Eigen::Array<float,3,Eigen::Dynamic> queryBatch(3,nPoints);
    Eigen::Array<float,4,Eigen::Dynamic> distBatch(4,nPoints);

    // Line of points perpindicular to the surface
    float ratio;
    for (int i = 0; i < nPoints; i++){
        ratio = (float)i/((float)nPoints-1.0);
        queryBatch(0,i) = 10.0;
        queryBatch(1,i) = 10.0;
        queryBatch(2,i) = -5.0 + 10.0*ratio;
    }

    // cout << "query batch is: " << queryBatch << endl;

    // Get distance
    distBatch = obj.getBatchDistanceFromPointsToSurface(queryBatch,101,101);

    // cout << "Distances to surface is: " << distBatch << endl;

    // Check output distances
    EXPECT_NEAR(5.0,distBatch(0,0),0.1);
    EXPECT_NEAR(-5.0,distBatch(0,nPoints-1),0.1);

    // Check output normals
    EXPECT_NEAR(0.0,distBatch(1,0),0.1);
    EXPECT_NEAR(0.0,distBatch(2,0),0.1);
    EXPECT_NEAR(-5.0,distBatch(3,0),0.1);
    
    EXPECT_NEAR(0.0,distBatch(1,nPoints-1),0.1);
    EXPECT_NEAR(0.0,distBatch(2,nPoints-1),0.1);
    EXPECT_NEAR(5.0,distBatch(3,nPoints-1),0.1);

}

TEST_F (object3DTest, getBatchDistanceToSurfaceAsymmetricObj){
    // Initialise class
    Object3D obj(scan2, 3, 3, 7, 10);

    // Test batch of points
    int nPoints = 35;
    Eigen::Array<float,3,Eigen::Dynamic> queryBatch(3,nPoints);
    Eigen::Array<float,4,Eigen::Dynamic> distBatch(4,nPoints);

    // Line of points perpindicular to the surface
    float ratio;
    for (int i = 0; i < nPoints; i++){
        ratio = (float)i/((float)nPoints-1.0);
        queryBatch(0,i) = 10.0;
        queryBatch(1,i) = 10.0;
        queryBatch(2,i) = -5.0 + 10.0*ratio;
    }

    // cout << "query batch is: " << queryBatch << endl;

    // Get distance
    distBatch = obj.getBatchDistanceFromPointsToSurface(queryBatch,101,101);

    // cout << "Distances to surface is: " << distBatch << endl;

    // Check output distances
    EXPECT_NEAR(5.0,distBatch(0,0),0.1);
    EXPECT_NEAR(-5.0,distBatch(0,nPoints-1),0.1);

    // Check output normals
    EXPECT_NEAR(0.0,distBatch(1,0),0.1);
    EXPECT_NEAR(0.0,distBatch(2,0),0.1);
    EXPECT_NEAR(-5.0,distBatch(3,0),0.1);
    
    EXPECT_NEAR(0.0,distBatch(1,nPoints-1),0.1);
    EXPECT_NEAR(0.0,distBatch(2,nPoints-1),0.1);
    EXPECT_NEAR(5.0,distBatch(3,nPoints-1),0.1);

}

class largeObject3DTest : public ::testing::Test {
 protected:
  largeObject3DTest() : scan(55,115) {
    float k = 0;
    
    for (int i = 0; i<scan.rows(); i++){
        for (int j=0; j<scan.cols(); j++){
            scan(i,j) = Point3Df(i,j,k);
        }
    }
  }

  ~largeObject3DTest() {;}


  Matrix_Point3Df scan;

};

TEST_F (largeObject3DTest, testGetSurfacePoints){
    Object3D obj(scan,3,3,30,30);

    // Number of points
    int msSurf = 25;
    int mtSurf = 25;

    // Init a point cloud
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>(mtSurf,msSurf,pcl::PointNormal())); 

    // Get the point cloud
    obj.getSurfacePointCloud(cloud,msSurf,mtSurf);

    float val(0.0001);

    // Test corner points
    EXPECT_FLOAT_EQ(cloud->at(0,0).x,obj.pointAt(val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,0).y,obj.pointAt(val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,0).z,obj.pointAt(val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).x,obj.pointAt(1.0-val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).y,obj.pointAt(1.0-val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).z,obj.pointAt(1.0-val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).x,obj.pointAt(val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).y,obj.pointAt(val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).z,obj.pointAt(val,1.0-val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).x,obj.pointAt(1.0-val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).y,obj.pointAt(1.0-val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).z,obj.pointAt(1.0-val,1.0-val).z());
}

TEST_F (largeObject3DTest, testGetSurfacePointsAsym){
    Object3D obj(scan,3,3,30,75);

    // Number of points
    int msSurf = 25;
    int mtSurf = 25;

    // Init a point cloud
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud (new pcl::PointCloud<pcl::PointNormal>(mtSurf,msSurf,pcl::PointNormal())); 

    // Get the point cloud
    obj.getSurfacePointCloud(cloud,msSurf,mtSurf);

    float val(0.0001);

    // Test corner points
    EXPECT_FLOAT_EQ(cloud->at(0,0).x,obj.pointAt(val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,0).y,obj.pointAt(val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,0).z,obj.pointAt(val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).x,obj.pointAt(1.0-val,val).x());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).y,obj.pointAt(1.0-val,val).y());
    EXPECT_FLOAT_EQ(cloud->at(0,msSurf-1).z,obj.pointAt(1.0-val,val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).x,obj.pointAt(val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).y,obj.pointAt(val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,0).z,obj.pointAt(val,1.0-val).z());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).x,obj.pointAt(1.0-val,1.0-val).x());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).y,obj.pointAt(1.0-val,1.0-val).y());
    EXPECT_FLOAT_EQ(cloud->at(mtSurf-1,msSurf-1).z,obj.pointAt(1.0-val,1.0-val).z());
}

TEST_F (largeObject3DTest, getDistanceToPoint){
    // Initialise class
    Object3D obj(scan,3,3,30,30);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    cout << "\n\nObj size is: (" << obj.ctrlPnts().rows() << ", " << obj.ctrlPnts().cols() << ")\n\n";

    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query,75,75);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(2.0,-dist,0.2);

}

TEST_F (largeObject3DTest, getDistanceToPointAsym){
    // Initialise class
    Object3D obj(scan,3,3,30,75);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    cout << "\n\nObj size is: (" << obj.ctrlPnts().rows() << ", " << obj.ctrlPnts().cols() << ")\n\n";


    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query, 75, 155);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(-2.0,dist,0.1);

}

TEST_F (largeObject3DTest, getDistanceToPointAsymlessTestPoints){
    // Initialise class
    Object3D obj(scan,3,3,30,75);

    // Declare point
    // Get distance to a point 
    Eigen::Vector3f query;    
    query << 10.0, 10.0, 2.0;
    cout << "query is " << query << endl;

    cout << "\n\nObj size is: (" << obj.ctrlPnts().rows() << ", " << obj.ctrlPnts().cols() << ")\n\n";


    // Get distance
    float dist = obj.getDistanceFromPointToSurface(query, 25, 25);

    cout << "Distance to surface is: " << dist << endl;

    EXPECT_NEAR(-2.0,dist,0.5);

}


} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
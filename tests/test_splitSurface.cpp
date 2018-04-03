
#include "cans/splitSurface.h"
#include "gtest/gtest.h"


namespace {

class getEndRowCountsTest : public ::testing::Test{
  protected:
    getEndRowCountsTest() : newPointsRows(5,1), newPointsCols(1,5), nRowColCount(1,4) {
        
        // COunt will be for L, R, D, U

        ss.nRowColCount.setZero(1,4);

    }
    ~getEndRowCountsTest(){;}

    Eigen::Array<int, 5, 1> newPointsRows;
    Eigen::Array<int, 1, 5> newPointsCols;
    Eigen::Array<int, 1, 4> nRowColCount;

    SplitSurface ss;

};

TEST_F (getEndRowCountsTest, testLeft){
    newPointsCols << 5, 5, 4, 0, 0;
    newPointsRows << 5, 0, 0, 0, 5;

    ss.newPointsRows = newPointsRows;
    ss.newPointsCols = newPointsCols;

    ss.getEndNewRowCounts();
    ss.getEndNewColCounts();

    // cout << "nRowColCount is [L, R, D, U]: [" << ss.nRowColCount << "]" << endl;    

    EXPECT_EQ(2, ss.nRowColCount[0]);
    EXPECT_EQ(0, ss.nRowColCount[1]);
    EXPECT_EQ(1, ss.nRowColCount[2]);
    EXPECT_EQ(1, ss.nRowColCount[3]);

}

TEST_F (getEndRowCountsTest, testRight){
    newPointsCols << 5, 0, 5, 5, 5;
    newPointsRows << 1, 0, 2, 0, 1;

    ss.newPointsRows = newPointsRows;
    ss.newPointsCols = newPointsCols;

    ss.getEndNewRowCounts();
    ss.getEndNewColCounts();

    // cout << "nRowColCount is [L, R, D, U]: [" << ss.nRowColCount << "]" << endl;    

    EXPECT_EQ(1, ss.nRowColCount[0]);
    EXPECT_EQ(3, ss.nRowColCount[1]);
    EXPECT_EQ(0, ss.nRowColCount[2]);
    EXPECT_EQ(0, ss.nRowColCount[3]);

}

TEST_F (getEndRowCountsTest, testDown){
    newPointsCols << 5, 0, 0, 3, 5;
    newPointsRows << 5, 5, 0, 1, 5;

    ss.newPointsRows = newPointsRows;
    ss.newPointsCols = newPointsCols;

    ss.getEndNewRowCounts();
    ss.getEndNewColCounts();

    // cout << "nRowColCount is [L, R, D, U]: [" << ss.nRowColCount << "]" << endl;    

    EXPECT_EQ(1, ss.nRowColCount[0]);
    EXPECT_EQ(1, ss.nRowColCount[1]);
    EXPECT_EQ(2, ss.nRowColCount[2]);
    EXPECT_EQ(1, ss.nRowColCount[3]);

}

TEST_F (getEndRowCountsTest, testUp){
    newPointsCols << 5, 0, 0, 5, 5;
    newPointsRows << 5, 1, 5, 5, 5;

    ss.newPointsRows = newPointsRows;
    ss.newPointsCols = newPointsCols;

    ss.getEndNewRowCounts();
    ss.getEndNewColCounts();

    // cout << "nRowColCount is [L, R, D, U]: [" << ss.nRowColCount << "]" << endl;    

    EXPECT_EQ(1, ss.nRowColCount[0]);
    EXPECT_EQ(2, ss.nRowColCount[1]);
    EXPECT_EQ(1, ss.nRowColCount[2]);
    EXPECT_EQ(3, ss.nRowColCount[3]);

}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

class splitSurfaceTest : public ::testing::Test {
 protected:
  splitSurfaceTest() : data(new pcl::PointCloud<pcl::PointNormal>(10,10)), data_rot(new pcl::PointCloud<pcl::PointNormal>(10,10)) 
  {
    
    int n_points = 10;

    // *data = pcl::PointCloud<pcl::PointNormal>(10,10);
    
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

    
    // Transform data to get test data sets   
    // Left 
    transform.translation() << -0.4, 0.2, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_l[0], transform);
    transform.translation() << -0.4, -0.2, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_l[1], transform);
    transform.translation() << -0.4, 0.0, 0.0;
    data_l.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_l[2], transform);  
    
    // Right
    transform.translation() << 0.4, 0.2, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_r[0], transform);
    transform.translation() << 0.4, -0.2, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_r[1], transform);
    transform.translation() << 0.4, 0.0, 0.0;
    data_r.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_r[2], transform);

    // Down
    transform.translation() << 0.2, -0.4, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_d[0], transform);
    transform.translation() << -0.2, -0.4, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_d[1], transform);
    transform.translation() << 0.0, -0.4, 0.0;
    data_d.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_d[2], transform);

    // Up
    transform.translation() << 0.2, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_u[0], transform);
    transform.translation() << -0.2, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_u[1], transform);
    transform.translation() << 0.0, 0.4, 0.0;
    data_u.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_u[2], transform);

    // Overlap (complete overlap)
    transform.scale(0.5f);// Scale by 1/2
    // cout << "transform for overlap is:\n" << transform.matrix() << endl;
    transform.translation() << 0.2, 0.2, 0.0;
    data_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_ol[0], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_ol[1], transform);
    transform.translation() << 0.4, 0.4, 0.0;
    data_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_ol[2], transform);

    // Mostly Overlap
    transform.translation() << -0.05, 0.2, 0.0;
    data_m_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_m_ol[0], transform);
    transform.translation() << 0.2, -0.05, 0.0;
    data_m_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_m_ol[1], transform);
    transform.translation() << 0.55, 0.55, 0.0;
    data_m_ol.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data, *data_m_ol[2], transform);  

    // Rotations
    Eigen::Vector3f rotVec(0.0, 0.0, 1); // Rotations about the z axis
    transform.scale(2.0f);// no scaling

    // Rotate 90 degrees (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (-M_PI/2.0f, rotVec));
    // cout << "transform for rot 90 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.4, 0.2, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[0], transform);
    transform.translation() << 0.4, 0.2, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[1], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[2], transform);
    transform.translation() << 0.2, -0.4, 0.0;
    data_rot1.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot1[3], transform);
    
    // Rotate 90 degrees   (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (M_PI, rotVec));// Rotate 180 from 90
    // cout << "transform for rot -90 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.4, 0.2, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[0], transform);
    transform.translation() << 0.4, 0.2, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[1], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[2], transform);
    transform.translation() << 0.2, -0.4, 0.0;
    data_rot2.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot2[3], transform);

    // Rotate 180 degrees   (-ve angle for a transformation)
    transform.rotate (Eigen::AngleAxisf (M_PI/2, rotVec));// Rotate another 90 to get to 180
    // cout << "transform for rot 180 is:\n" << transform.matrix() << endl;
    transform.translation() << -0.4, 0.2, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot3[0], transform);
    transform.translation() << 0.4, 0.2, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot3[1], transform);
    transform.translation() << 0.2, 0.4, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
    pcl::transformPointCloud(*data_rot, *data_rot3[2], transform);
    transform.translation() << 0.2, -0.4, 0.0;
    data_rot3.push_back(pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>(10,10)));
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
  // For rotation
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot1;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot2;
  std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr, Eigen::aligned_allocator <pcl::PointCloud <pcl::PointNormal>::Ptr > > data_rot3;  

  // Class instatiation
  SplitSurface ss;

};

// Constructors

TEST_F (splitSurfaceTest, BlankConstructor){
    SplitSurface ss;
}

TEST_F (splitSurfaceTest, InputMapCloud){
    ss.setInputMap(data);
}

TEST_F (splitSurfaceTest, InputObservationCloud){
    ss.setInputObservation(data);
}

TEST_F (splitSurfaceTest, InputObservationCloudVec){
    ss.setInputObservation(data_l[0]);
}

TEST_F (splitSurfaceTest, distThresholdCalcs){

    ss.setInputMap(data);

    float dt = ss.computeDistanceThreshold();

    cout << "Distance theshold is 0.3*" << dt << endl;
}

// Left test
TEST_F (splitSurfaceTest, leftTests){

  std::string expectExtDir = "LL";

  Eigen::Array<int, 2, 2> expectNewIndices;
  expectNewIndices << 0, 0,
                      9, 2; 

  Eigen::Array<int, 2, 2> expectOverlapIndices;
  expectOverlapIndices << 0, 2,
                          9, 9; 



  for (int i = 0; i < 3; i++){
    // Run split surface
    ss.splitNewSurfaceObservation(data,data_l[i]);


    // Check result
    EXPECT_EQ(expectExtDir, ss.extendDirection);
    EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
  }

}

// Right test
TEST_F (splitSurfaceTest, rightTests){

  std::string expectExtDir = "RR";

  Eigen::Array<int, 2, 2> expectNewIndices;
  expectNewIndices << 0, 7,
                      9, 9; 

  Eigen::Array<int, 2, 2> expectOverlapIndices;
  expectOverlapIndices << 0, 0,
                          9, 7; 



  for (int i = 0; i < 3; i++){
    // Run split surface
    ss.splitNewSurfaceObservation(data,data_r[i]);


    // Check result
    EXPECT_EQ(expectExtDir, ss.extendDirection);
    EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
  }

}

// Down test
TEST_F (splitSurfaceTest, downTests){

  std::string expectExtDir = "DD";

  Eigen::Array<int, 2, 2> expectNewIndices;
  expectNewIndices << 0, 0,
                      2, 9; 

  Eigen::Array<int, 2, 2> expectOverlapIndices;
  expectOverlapIndices << 2, 0,
                          9, 9; 



  for (int i = 0; i < 3; i++){
    // Run split surface
    ss.splitNewSurfaceObservation(data,data_d[i]);


    // Check result
    EXPECT_EQ(expectExtDir, ss.extendDirection);
    EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
  }

}

// Up test
TEST_F (splitSurfaceTest, upTests){

  std::string expectExtDir = "UU";

  Eigen::Array<int, 2, 2> expectNewIndices;
  expectNewIndices << 7, 0,
                      9, 9; 

  Eigen::Array<int, 2, 2> expectOverlapIndices;
  expectOverlapIndices << 0, 0,
                          7, 9; 

  for (int i = 0; i < 3; i++){
    // Run split surface
    ss.splitNewSurfaceObservation(data,data_u[i]);

    // Check result
    EXPECT_EQ(expectExtDir, ss.extendDirection);
    EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
  }
}

TEST_F (splitSurfaceTest, olTest){

  std::string expectExtDir = "NN";

  Eigen::Array<int, 2, 2> expectNewIndices;
  expectNewIndices << 0, 0,
                      0, 0; 

  Eigen::Array<int, 2, 2> expectOverlapIndices;
  expectOverlapIndices << 0, 0,
                          9, 9; 

  for (int i = 0; i < 3; i++){
    // Run split surface
    ss.splitNewSurfaceObservation(data,data_ol[i]);

    // Check result
    EXPECT_EQ(expectExtDir, ss.extendDirection);
    EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
  }    
}

TEST_F (splitSurfaceTest, mostlyOlTest){

  std::string expectExtDir = "NN";

  Eigen::Array<int, 2, 2> expectNewIndices;
  expectNewIndices << 0, 0,
                      0, 0; 

  Eigen::Array<int, 2, 2> expectOverlapIndices;
  expectOverlapIndices << 0, 0,
                          9, 9; 

  for (int i = 0; i < 3; i++){
    // Run split surface
    ss.splitNewSurfaceObservation(data,data_m_ol[i]);

    // Check result
    EXPECT_EQ(expectExtDir, ss.extendDirection);
    EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
  }    
}

TEST_F (splitSurfaceTest, rotate90Test){
    
  // initialise
  std::string expectExtDir = "NN"; 
  Eigen::Array<int, 2, 2> expectNewIndices;
  Eigen::Array<int, 2, 2> expectOverlapIndices;

  for (int i = 0; i < 4; i ++){
    // Run split surface
    ss.splitNewSurfaceObservation(data_rot,data_rot1[i]);

    switch (i) {
      case 0 : 
        expectExtDir = "LD";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 0,
                            2, 9; 
        expectOverlapIndices << 2, 0,
                                9, 9;                     
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 1 : 
        expectExtDir = "RU";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 7, 0,
                            9, 9; 
        expectOverlapIndices << 0, 0,
                                7, 9;                     
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 2 : 
        expectExtDir = "UL";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 0,
                            9, 2; 
        expectOverlapIndices << 0, 2,
                                9, 9;                     
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 3 : 
        expectExtDir = "DR";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 7,
                            9, 9; 
        expectOverlapIndices << 0, 0,
                                9, 7;                     
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
    }
  }
}

TEST_F (splitSurfaceTest, rotateMinus90Test){
    
  // initialise
  std::string expectExtDir = "NN"; 
  Eigen::Array<int, 2, 2> expectNewIndices;
  Eigen::Array<int, 2, 2> expectOverlapIndices;

  for (int i = 0; i < 4; i ++){
    // Run split surface
    ss.splitNewSurfaceObservation(data_rot,data_rot2[i]);

    switch (i) {
      case 0 : 
        expectExtDir = "LU";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 7, 0,
                            9, 9; 
        expectOverlapIndices << 0, 0,
                                7, 9;         
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 1 : 
        expectExtDir = "RD";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 0,
                            2, 9; 
        expectOverlapIndices << 2, 0,
                                9, 9;                                         
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 2 : 
        expectExtDir = "UR";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 7,
                            9, 9; 
        expectOverlapIndices << 0, 0,
                                9, 7;
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 3 : 
        expectExtDir = "DL";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 0,
                            9, 2; 
        expectOverlapIndices << 0, 2,
                                9, 9;  
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
    }
  }
}

TEST_F (splitSurfaceTest, rotate180Test){
    
  // initialise
  std::string expectExtDir = "NN"; 
  Eigen::Array<int, 2, 2> expectNewIndices;
  Eigen::Array<int, 2, 2> expectOverlapIndices;

  for (int i = 0; i < 4; i ++){
    // Run split surface
    ss.splitNewSurfaceObservation(data_rot,data_rot3[i]);

    switch (i) {
      case 0 : 
        expectExtDir = "LR";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 7,
                            9, 9; 
        expectOverlapIndices << 0, 0,
                                9, 7;
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 1 : 
        expectExtDir = "RL";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 0,
                            9, 2; 
        expectOverlapIndices << 0, 2,
                                9, 9;  
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 2 : 
        expectExtDir = "UD";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 0, 0,
                            2, 9; 
        expectOverlapIndices << 2, 0,
                                9, 9;
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
      case 3 : 
        expectExtDir = "DU";
        EXPECT_EQ(expectExtDir, ss.extendDirection);
        expectNewIndices << 7, 0,
                            9, 9; 
        expectOverlapIndices << 0, 0,
                                7, 9;         
        EXPECT_TRUE((expectNewIndices == ss.newDataIndices).all());
        break;
    }
  }
}

} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
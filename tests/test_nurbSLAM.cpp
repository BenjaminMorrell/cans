#include "cans/nurbSLAM.h"
#include "gtest/gtest.h"

namespace {

class nurbSLAMIOTest : public ::testing::Test{
  protected:
    nurbSLAMIOTest() {
      ;
    }

    ~nurbSLAMIOTest() {;}


};

TEST_F (nurbSLAMIOTest, constructor){
  NurbSLAM ns;
}

TEST_F (nurbSLAMIOTest, loadObject){
  NurbSLAM ns;

  std::string filename = "/home/bjm/SpaceCRAFT/ros_ws/src/cans/examples/blob_scan_res_new_extend_final_obj.obj";

  ns.loadObjectIntoMap(filename);

  EXPECT_EQ(1,ns.mp.objectMap.size());
  // EXPECT_EQ(1,ns.mapMeshList.size());
  // EXPECT_EQ(1,ns.mapFeatureList.size());
}

TEST_F (nurbSLAMIOTest, loadTwoObjects){
  NurbSLAM ns;

  std::string filename = "/home/bjm/SpaceCRAFT/ros_ws/src/cans/examples/blob_scan_res_new_extend_final_obj.obj";

  ns.loadObjectIntoMap(filename);

  ns.loadObjectIntoMap(filename);

  EXPECT_EQ(2,ns.mp.objectMap.size());
  // EXPECT_EQ(2,ns.mapMeshList.size());
  // EXPECT_EQ(2,ns.mapFeatureList.size());
}


// class nurbSLAMTest : public ::testing::Test{
//   protected:
//     nurbSLAMTest() {
//       ;
//     }

//     ~nurbSLAMTest() {;}

//     NurbSLAM ns;

// };

} // namespace

int main (int argc,char ** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#ifndef _cans_nurbSLAM_h_
#define _cans_nurbSLAM_h_
// This guards the headers so you don't accidently include them twice

#include "cans/mapping3D.h"

//PCL includes
#include <Eigen/Core>
// Already in mapping3d
// #include <pcl/io/pcd_io.h>
// #include <pcl/point_types.h>
// #include <pcl/point_cloud.h>
// #include <pcl/common/time.h>
// #include <pcl/common/transforms.h>

// Already in splitSurface.h
// #include <pcl/features/normal_3d_omp.h>
// #include <pcl/features/fpfh_omp.h>
#include <pcl/features/normal_3d.h>

// Already in splitSurface.h
// #include <pcl/registration/icp.h>
// #include <pcl/registration/sample_consensus_prerejective.h>
// #include <pcl/registration/correspondence_estimation_normal_shooting.h>
// #include <pcl/registration/correspondence_rejection_distance.h>
#include <pcl/registration/ia_ransac.h>
#include <pcl/segmentation/sac_segmentation.h>
// May not need..
#include <pcl/visualization/pcl_visualizer.h>
typedef pcl::visualization::PointCloudColorHandlerCustom<pcl::PointNormal> ColorHandlerT; // Visualisation?


class NurbSLAM {

  private:
    
    Eigen::Affine3f state;
    Eigen::Affine3f transformDelta;

  
    // Other SLAM filter stuff maybe...

    // Lists to store sets of scans
    std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> objectMeshList; // This may break
    std::vector<int> objIDList;
    std::vector<Eigen::Matrix4f> transformationList;

  public:

    NurbSLAM();

    ~NurbSLAM();

    // Higher level function
    void processScans(std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds);

    int processSingleScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed);

    Eigen::Matrix4f alignScanWithMapObject(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloud2);

    void updateSLAMFilter();

    void alignAndUpdateMeshes();

    void initState(Eigen::Affine3f startingState);

    Eigen::Affine3f getState();

    // Map + Mapping functions
    Mapping3D mp;

    bool bShowAlignment;

    // Localisation settings
    float nSurfPointsFactor;// - factor multiplied by the number of control points to get the number of surface samples
    float pclRadiusSetting; // size to search for normal and feature estimation.

    int localisationOption;
};


#endif
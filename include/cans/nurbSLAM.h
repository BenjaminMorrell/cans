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

#include <pcl/keypoints/iss_3d.h>
#include <pcl/keypoints/harris_3d.h>
#include <pcl/keypoints/harris_6d.h>
#include <pcl/keypoints/smoothed_surfaces_keypoint.h>

// May not need..
#include <pcl/visualization/pcl_visualizer.h>
typedef pcl::visualization::PointCloudColorHandlerCustom<pcl::PointNormal> ColorHandlerT; // Visualisation?


class NurbSLAM {

  private:
    
    Eigen::Affine3f state;
    Eigen::Affine3f transformDelta;

  
    // Other SLAM filter stuff maybe...

    // Lists to store sets of scans
    std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> objectMeshList; 
    std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> mapMeshList; 
    std::vector<pcl::PointCloud<pcl::FPFHSignature33>::Ptr> mapFeatureList; 
    std::vector<int> objIDList;
    std::vector<Eigen::Matrix4f> transformationList;

    float inlierFraction;
    

    // Options
    bool bMappingModeOn;
    bool bLocalisationModeOn;

  public:

    NurbSLAM();

    ~NurbSLAM();

    // Higher level function
    void processScans(std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds);

    int processSingleScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed);

    Eigen::Matrix4f alignScanKeypointsWithMapObjectKeypoints(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);
    Eigen::Matrix4f alignScanKeypointsWithMapObjectDense(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);
    Eigen::Matrix4f alignScanWithMapObject(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);
    void computeKeypoints(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_, pcl::PointCloud<pcl::PointNormal>::Ptr keypoints);

    void updateSLAMFilter();

    void alignAndUpdateMeshes();
    void updatePointCloudAndFeaturesInMap(int objID);

    void setState(Eigen::Affine3f startingState);
    Eigen::Affine3f getState();

    // Changing modes
    void activateMappingMode();
    void activateLocalisationMode();

    bool isMappingModeActive(){return bMappingModeOn;}
    bool isLocalisationModeActive(){return bLocalisationModeOn;}

    // Load in map
    void loadObjectIntoMap(std::string filename);

    // Map + Mapping functions
    Mapping3D mp;

    // Options
    bool bShowAlignment;    

    // Flags
    bool bMapUpdatedFromScan; // Flag for ros node to know if the map has changed

    // Localisation settings
    float nSurfPointsFactor;// - factor multiplied by the number of control points to get the number of surface samples
    float pclNormalRadiusSetting; // size to search for normal estimation.
    float pclFeatureRadiusSetting; // size to search for feature estimation.
    float validInlierTheshold;
    
    int alignmentOption;
    int localisationOption;

    float modelResolutionKeypoints; // Setting for keypoint extraction
    int minNeighboursKeypoints;

    float ransac_inlierMultiplier; // Setting that affects the size of the inlier threshold
    int ransac_maximumIterations;
    int ransac_numberOfSamples;
    int ransac_correspondenceRandomness;
    float ransac_similarityThreshold;
    float ransac_inlierFraction;


    int keypointOption;
};


#endif
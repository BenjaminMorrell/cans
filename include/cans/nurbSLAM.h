#ifndef _cans_nurbSLAM_h_
#define _cans_nurbSLAM_h_
// This guards the headers so you don't accidently include them twice

#include "cans/mapping3D.h"
#include <chrono>

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
#include <pcl/filters/extract_indices.h>
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
    Eigen::Affine3f oldState; // after the process step
    Eigen::Affine3f previousState; // Before most recent process step
    Eigen::Affine3f transformDelta;

  
    // SLAM components
    Eigen::Matrix<float,12,1> ekfState;
    Eigen::Matrix<float,12,12> P;
    Eigen::Matrix<float,12,12> J;
    Eigen::Matrix<float,12,12> Q;
    Eigen::Matrix<float,6,6> R;
    Eigen::Matrix<float,6,12> Jh;
    Eigen::Matrix<float,12,6> K;

    // Lists to store sets of scans
    std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> objectMeshList; 
    std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> mapMeshList; 
    std::vector<pcl::PointCloud<pcl::FPFHSignature33>::Ptr> mapFeatureList; 
    std::vector<int> objIDList;
    std::vector<Eigen::Matrix4f> transformationList;
    std::vector<float> inlierFractionList;
    std::vector<int> mapMatchCount;
    std::vector<int> mapExtendCount;

    float inlierFraction;
    int numberOfPointsInAlignment;
    bool bRejectAlignment;
    
    bool bObjectNormalsComputed;

    // Options
    bool bMappingModeOn;
    bool bLocalisationModeOn;

    Object3D oldObj;

  public:

    NurbSLAM();

    ~NurbSLAM();

    // Higher level function
    void processScans(std::vector<pcl::PointCloud<pcl::PointNormal>::Ptr> clouds, float timestep = 0.0);

    int processSingleScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr cloudTransformed);

    Eigen::Matrix4f alignScanKeypointsWithMapObjectKeypoints(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);
    Eigen::Matrix4f alignScanKeypointsWithMapObjectDense(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);
    Eigen::Matrix4f alignScanWithMapObject(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);
    void computeKeypoints(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::search::KdTree<pcl::PointNormal>::Ptr search_method_, pcl::PointCloud<pcl::PointNormal>::Ptr keypoints);
    void computeNormals(pcl::PointCloud<pcl::PointNormal>::Ptr cloud);
    void computeFeatures(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, pcl::PointCloud<pcl::PointNormal>::Ptr searchSurface, pcl::PointCloud<pcl::FPFHSignature33>::Ptr features);
    void rejectNonOverlappingPoints(pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsPCFilt);

    void processStepEKF(float timestep);
    void updateSLAMFilter(float timestep);

    void alignAndUpdateMeshes();
    void updatePointCloudAndFeaturesInMap(int objID);

    void setInitEKFStates();

    void setState(Eigen::Affine3f startingState);
    Eigen::Affine3f getState();

    bool doesPointCloudHaveNans(pcl::PointCloud<pcl::PointNormal>::Ptr cloud);

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

    bool bRejectNonOverlappingInAlign; // To reject points that don't overlap before alignment
    float maxDistanceOverlap;
    
    int mapCountThreshold;
    int mapExtendThreshold;

    int alignmentOption;
    int localisationOption;

    bool bUseFullAlignmentTransformInUpdate;
    bool bUseOldStateForNewObjects; 
    std::vector<float> rejectCriteria;
    bool bKeepPConstant;

    float modelResolutionKeypoints; // Setting for keypoint extraction
    int minNeighboursKeypoints;

    float ransac_inlierMultiplier; // Setting that affects the size of the inlier threshold
    int ransac_maximumIterations;
    int ransac_numberOfSamples;
    int ransac_correspondenceRandomness;
    float ransac_similarityThreshold;
    float ransac_inlierFraction;

    std::vector<double> processTimes; // Vector that is 5 values long. 1) Mesh processing, 2) Data association, 3) Alignment, 4) pose update, 5) Update map
    // For a given process Scan call, this is cumulative for the number of scans

    int keypointOption;

    // SLAM options

    float pNoisePos;
    float pNoiseVel;
    float pNoiseAccel;
    float pNoiseAng;
    float pNoiseMultiplier;
    float qNoiseMultiplier;
    float noiseObsBasePos;
    float noiseObsMultPos;
    float noiseObsMultPosErr;
    float noiseObsBaseAng;
    float noiseObsMultAng;
    float noiseObsMultAngErr;
    float rMatMultiplier;

    int processModel; // 0 - const accel, 1 - const vel, 2 - const pos
};


#endif

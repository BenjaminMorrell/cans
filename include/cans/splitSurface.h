#ifndef _cans_SplitSurface_h_
#define _cans_SplitSurface_h_

#include <nurbsS.h>
// #include <cmath>
// PCL specific includes
// #include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/fpfh_omp.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/sample_consensus_prerejective.h>
#include <pcl/registration/correspondence_estimation_normal_shooting.h>
#include <pcl/registration/correspondence_estimation.h>
#include <pcl/registration/correspondence_rejection_distance.h>


class SplitSurface {

  private:

    pcl::CorrespondencesPtr correspondences_;
    // pcl::CorrespondencesPtr corr_filtPtr;
    
    float maxDistThreshMultiplier;

  public:

    SplitSurface();

    ~SplitSurface();

    void splitNewSurfaceObservation(pcl::PointCloud<pcl::PointNormal>::Ptr mapCloud, pcl::PointCloud<pcl::PointNormal>::Ptr obsCloud);

    void setInputMap(pcl::PointCloud<pcl::PointNormal>::Ptr);

    void setInputObservation(pcl::PointCloud<pcl::PointNormal>::Ptr);

    void findNonOverlappingData();

    float computeDistanceThreshold();

    void getNewDataExtendDirection ();

    void getEndNewRowCounts();

    void getEndNewColCounts();

    void getMapDataExtendDirection();

    int getCloudWidth();

    

    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> newPointsArray;
    std::string extendDirection = "NN";

    Eigen::Array<int, 2,2> newDataIndices;
    Eigen::Array<int, 2, 2> overlapDataIndices;

    Eigen::Array<int, Eigen::Dynamic, 1> newPointsRows;
    Eigen::Array<int, 1, Eigen::Dynamic> newPointsCols;

    Eigen::Array<int, 1, 4> nRowColCount;

    int nExtraNew; // Number of new rows or columns to take
    int newRowColBuffer; // The allowed number of overlap points in a row or column

    pcl::PointCloud<pcl::PointNormal>::Ptr mapCloud;
    pcl::PointCloud<pcl::PointNormal>::Ptr obsCloud;

};


#endif
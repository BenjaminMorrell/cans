
#ifndef _cans_Mapping3D_h_
#define _cans_Mapping3D_h_

#include <nurbsS.h>
// #include <cmath>
// PCL specific includes
// #include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

using namespace PLib ; 

class Mapping3D {

private:
    int n_ctrl_default;
    int order[2]; // Work out how to do int vectors...
    int number_downsample_points;

    float overlap_d_frac;

    float search_thresh[9];




public:
    Mapping3D() ;
    // Other constructors here?
    ~Mapping3D(){;}

    // std::vector<float> compute_centre_of_data()
    // Or use Eigen - for handling matrices
    // Eigen::MatrixXd 
    // PLib uses Eigen anyway

    NurbsCurvef joinCurves(NurbsCurvef& crv1, NurbsCurvef& crv2, bool newBeforeOld = false, bool flipKnotParam = false);
    NurbsSurfacef joinSurfaces(NurbsSurfacef&, NurbsSurfacef&, char *);

    // Scan processing functions
    void meshFromScan(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOut, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudIn);
    void getNanMatrixFromPointCloud(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);
    bool removeRowsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags);
    bool removeColsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags);
    void downsampleRow(Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags);
    void downsampleCol(Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags);
    void averageOutNans(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic> nanIndices);
    void regionAverage(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, int i, int j);
    
    // Convenience functions
    Vector_HPoint3Df getMatRow(Matrix_HPoint3Df, int);
    Vector_HPoint3Df getMatCol(Matrix_HPoint3Df, int, bool);
    void insertMatRow(Matrix_HPoint3Df&,Vector_HPoint3Df, int);
    void insertMatCol(Matrix_HPoint3Df&,Vector_HPoint3Df, int, bool);


    int numRowsDesired;
    int numColsDesired;
    int maxNanAllowed;
    int removeNanBuffer;


};

#endif
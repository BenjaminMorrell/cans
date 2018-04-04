
#ifndef _cans_Mapping3D_h_
#define _cans_Mapping3D_h_

#include "cans/object3D.h"
#include "cans/splitSurface.h"
#include <nurbsS.h>
#include <cmath>
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

    std::vector<float> searchThresh;

    
    int numberOfMetrics;



public:
    Mapping3D() ;
    // Other constructors here?
    ~Mapping3D(){;}

    // std::vector<float> compute_centre_of_data()
    // Or use Eigen - for handling matrices
    // Eigen::MatrixXd 
    // PLib uses Eigen anyway

    // Surface extension functions
    NurbsCurvef joinCurves(NurbsCurvef& crv1, NurbsCurvef& crv2, bool newBeforeOld = false, bool flipKnotParam = false);
    NurbsSurfacef joinSurfaces(NurbsSurfacef&, NurbsSurfacef&, std::string);
    Object3D joinSurfaces(Object3D& srf1, Object3D& srf2, std::string extendDir);
    std::vector<int> computeNumberOfControlPoints(std::string extendDirection, Matrix_Point3Df& data, Object3D& srf);
    
    // Scan processing functions
    void meshFromScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloudOut, pcl::PointCloud<pcl::PointNormal>::Ptr cloudIn);
    void getNanMatrixFromPointCloud(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray, pcl::PointCloud<pcl::PointNormal>::Ptr cloud);
    bool removeRowsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags);
    bool removeColsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags);
    void downsampleRow(Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags);
    void downsampleCol(Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags);
    bool averageOutNans(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic>& nanIndices);
    void regionAverage(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, int i, int j);

    // Object Update
    void addObject(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, std::vector<float> searchMetrics);
    void updateObject(Object3D& mapObj, pcl::PointCloud<pcl::PointNormal>::Ptr mapObjPC, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);


    // Data association
    int dataAssociation(std::vector<float> searchMetrics);

    
    // Convenience functions
    Matrix_Point3Df nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);
    Matrix_Point3Df nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud);
    Matrix_Point3Df nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int, 2, 2> dataIndices);
    pcl::PointCloud<pcl::PointNormal> pointCloudFromNurbsData(Matrix_Point3Df data);
    Vector_HPoint3Df getMatRow(Matrix_HPoint3Df, int);
    Vector_HPoint3Df getMatCol(Matrix_HPoint3Df, int, bool);
    void insertMatRow(Matrix_HPoint3Df&,Vector_HPoint3Df, int);
    void insertMatCol(Matrix_HPoint3Df&,Vector_HPoint3Df, int, bool);


    int numRowsDesired;
    int numColsDesired;
    int maxNanAllowed;
    int removeNanBuffer;

    std::vector<Object3D> objectMap;
    // Eigen::Array<float,Eigen::Dynamic, Eigen::Dynamic> objectMetrics;
    // or 
    std::vector<std::vector<float> > objectMetrics;


};

#endif
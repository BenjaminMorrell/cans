
#ifndef _cans_Mapping3D_h_
#define _cans_Mapping3D_h_

#include "cans/object3D.h"
#include "cans/splitSurface.h"
#include <nurbsS.h>

#include <cans_msgs/Object3D.h>

// PCL specific includes
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/common/transforms.h>
#include <pcl/common/spring.h>



using namespace PLib ; 

class Mapping3D {

private:
    
    

    

    bool knotInsertionFlag;
    int numInsert;
    float deltaKnotInsert;
    
    int numberOfMetrics;


public:
    Mapping3D() ;
    // Other constructors here?
    ~Mapping3D(){;}

    // std::vector<float> compute_centre_of_data()
    // Or use Eigen - for handling matrices
    // Eigen::MatrixXd 
    // PLib uses Eigen anyway

    // Higher level operations
    int processScan(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Affine3f transform);

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
    void downsampleRow(Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags, int numRows = -1);
    void downsampleCol(Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags, int numCols = -1);
    bool averageOutNans(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic>& nanIndices);
    void regionAverage(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, int i, int j);
    void regionAverageSimple(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, int i, int j);

    // Object Update
    void addObject(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, std::vector<float> searchMetrics);
    void addObjectFromFile(const char * filename);
    void updateObjectInMap(int objID, Object3D& obj);
    void updateObject(int objID, pcl::PointCloud<pcl::PointNormal>::Ptr obsObjPC);
    void knotInsertionPreMerge(Object3D& obj, std::string extendDirection);
    void knotInsertionAlongSeam(Object3D& obj, std::string extendDirection, int nInsert);

    // Data association
    std::vector<float> computeSearchMetrics(Object3D& obj);
    std::vector<float> computeSearchMetrics(pcl::PointCloud<pcl::PointNormal>::Ptr cloud);
    int dataAssociation(std::vector<float> searchMetrics);

    // Convenience functions
    Matrix_Point3Df nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);
    Matrix_Point3Df nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud);
    Matrix_Point3Df nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int, 2, 2>& dataIndices);
    Matrix_Point3Df nurbsDataFromPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr cloud, Eigen::Array<int, Eigen::Dynamic, 2>& newRowIndices, Eigen::Array<int, Eigen::Dynamic, 2>& newColIndices);
    void pointCloudFromNurbsData(Matrix_Point3Df& data, pcl::PointCloud<pcl::PointNormal>::Ptr cloud);
    void pointCloudFromObject3D(int objID, int ms, int mt, pcl::PointCloud<pcl::PointNormal>::Ptr cloud);
    Vector_HPoint3Df getMatRow(Matrix_HPoint3Df, int);
    Vector_HPoint3Df getMatCol(Matrix_HPoint3Df, int, bool);
    void insertMatRow(Matrix_HPoint3Df&,Vector_HPoint3Df, int);
    void insertMatCol(Matrix_HPoint3Df&,Vector_HPoint3Df, int, bool);

    int getNumberOfUniquePoints(Eigen::Array<int, Eigen::Dynamic, 2>& Indices);

    // Output functions
    void writeObjectPCDFile(const char* filename, const int objID, int ms = -1, int mt = -1);

    void fillObject3DMessage(int objID, cans_msgs::Object3D& msg);

    SplitSurface ss;

    int numRowsDesired;
    int numColsDesired;
    int minRowsColsAllowed;
    int maxNanAllowed;
    float maxNanPercentage;
    int removeNanBuffer;
    int nCtrlDefault[2];
    int order[2]; 

    std::vector<float> searchThresh;

    bool useNonRectData;

    bool bFilterZ;
    int nPointsZLim;
    float zThreshMultiplier;
    bool bRejectScan;
    
    bool bNegateZ;

    int msSurf;
    int mtSurf;

    int newRowColBuffer; // The allowed number of overlap points in a row or column

    std::vector<Object3D> objectMap;
    // Eigen::Array<float,Eigen::Dynamic, Eigen::Dynamic> objectMetrics;
    // or 
    std::vector<std::vector<float> > objectMetrics;


};

#endif
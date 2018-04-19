
#ifndef _cans_object3D_h_
#define _cans_object3D_h_
// This guards the headers so you don't accidently include them twice

#include <nurbsS.h>
#include <nurbs.h>

#include <Eigen/Core>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/registration/correspondence_estimation.h>

namespace PLib{

class Object3D : public NurbsSurfacef {
  protected:
    Point3Df centre;
    Point3Df color;
    float objSize;
    float GISTexture;
    float GISRoughness;

  public:
    // Empty, default
    Object3D() ;

    // Copy constructor 
    Object3D(const Object3D & obj);

    // Surface input with defaults for properties
    // Object3D(PlNurbsSurfacef &nS, Point3Df &cp, Point3Df &col, float objSize = 0, float GIS1 = 0.0, float GIS2 = 0.0) ;

    // LS fit constructor
    Object3D(const Matrix_Point3Df& Q); // with defaults
    Object3D(const Matrix_Point3Df& Q, int pU, int pV, int nU, int nV);

    // Knots and control points constructor
    Object3D(int pU, int pV, Vector_FLOAT& Uvec, Vector_FLOAT& Vvec, Matrix_HPoint3Df& ctrlPnts);
    // using Eigen inputs
    Object3D(int pU, int pV, Eigen::Array<float,1,Eigen::Dynamic>& Uvec, Eigen::Array<float,1,Eigen::Dynamic>& Vvec, Eigen::Array<float,3,Eigen::Dynamic>& ctrlPnts, int nCtrlS, int nCtrlT);

    ~Object3D();

    void operator = (const Object3D& obj);
    
    // virtual ~Object3D(){;}// empty destructor?
  public:

    void updateObject3D(int pU, int pV, Eigen::Array<float,1,Eigen::Dynamic>& Uvec, Eigen::Array<float,1,Eigen::Dynamic>& Vvec, Eigen::Array<float,3,Eigen::Dynamic>& ctrlPnts, int nCtrlS, int nCtrlT);


    void computeCentreFromData(const Matrix_Point3Df& scan, int step_size = 1);
    void computeCentreFromControlPoints();


    void computeSizeFromData(const Matrix_Point3Df& Q);
    void computeSizeFromControlPoints();

    void readObject3D(const char* filename);

    Matrix_Point3Df getSurfacePoints(int ms = 45, int mt = 45);
    void getSurfacePointCloud( pcl::PointCloud<pcl::PointNormal>::Ptr cloud,int ms = 45, int mt = 45);

    float getDistanceFromPointToSurface(Eigen::Vector3f& query, int ms = 20, int mt = 20);
    Eigen::Array<float,1,Eigen::Dynamic> getBatchDistanceFromPointsToSurface(Eigen::Array<float,3,Eigen::Dynamic>& query, int ms = 20, int mt = 20);



    Point3Df& getCentre();
    Point3Df& getColor();
    float getObjSize();
    // float getGISTexture();
    // float getGISRoughness();

};


} // Namespace


#endif 

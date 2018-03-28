
#ifndef _cans_object3D_h_
#define _cans_object3D_h_
// This guards the headers so you don't accidently include them twice

#include <nurbsS.h>
#include <nurbs.h>

namespace PLib{

class Object3D : public NurbsSurfacef {
  protected:
    Point3Df centre;
    Point3Df color;
    float obj_size;
    float GIS_texture;
    float GIS_roughness;

  public:
    // Empty, default
    Object3D() ;

    // Surface input with defaults for properties
    // Object3D(PlNurbsSurfacef &nS, Point3Df &cp, Point3Df &col, float obj_size = 0, float GIS1 = 0.0, float GIS2 = 0.0) ;

    // LS fit constructor
    Object3D(const Matrix_Point3Df& Q); // with defaults
    Object3D(const Matrix_Point3Df& Q, int pU, int pV, int nU, int nV);
    
    // virtual ~Object3D(){;}// empty destructor?
  public:

    void computeCentreFromData(const Matrix_Point3Df& scan, int step_size = 1);
    void computeCentreFromControlPoints();


    void computeSizeFromData(const Matrix_Point3Df& Q);
    void computeSizeFromControlPoints();

    Point3Df& getCentre();
    Point3Df& getColor();
    float getObjSize();
    // float getGISTexture();
    // float getGISRoughness();

};


} // Namespace


#endif 

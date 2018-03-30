
#ifndef _cans_Mapping3D_h_
#define _cans_Mapping3D_h_

#include <nurbsS.h>

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

    Point3Df compute_centre_of_data(const Matrix_Point3Df&);

    // std::vector<float> compute_centre_of_data()
    // Or use Eigen - for handling matrices
    // Eigen::MatrixXd 
    // PLib uses Eigen anyway

    Matrix_Point3Df mesh_from_scan(const Matrix_Point3Df);

    NurbsCurvef joinCurves(NurbsCurvef& crv1, NurbsCurvef& crv2, bool newBeforeOld = false, bool flipKnotParam = false);
    NurbsSurfacef joinSurfaces(NurbsSurfacef&, NurbsSurfacef&, char *);

    // Convenience functions
    Vector_HPoint3Df getMatRow(Matrix_HPoint3Df, int);
    Vector_HPoint3Df getMatCol(Matrix_HPoint3Df, int, bool);
    void insertMatRow(Matrix_HPoint3Df&,Vector_HPoint3Df, int);
    void insertMatCol(Matrix_HPoint3Df&,Vector_HPoint3Df, int, bool);

};

#endif
#include "ContourMethod.h"
#include "Inflection_Curve.h"

int Inflection_Curve::function_on_square(double *foncub, int i, int j) {
    int is_square = gv->cell_type(i, j);

    double rvorig[2], rvside[2], rvtop[2], rvopps[2];
    double forig, fside, ftop, fopps;
    int orient;

    // Test for real eigenvalues on the triangle
    if (!gv->eig_is_real(i, j)[family]) return 0;
    if (!gv->eig_is_real(i + 1, j)[family]) return 0;
    if (!gv->eig_is_real(i, j + 1)[family]) return 0;

    for (int k = 0; k < 2; k++) {
        rvorig[k] = gv->e(i, j)[family].vrr[k];
        rvside[k] = gv->e(i + 1, j)[family].vrr[k];
        rvtop[k]  = gv->e(i, j + 1)[family].vrr[k];
    }

    // Make vectors consistent at pairs of corners (if possible)
    // if vectors consistent, compute directional derivatives
    if (consistency(rvtop, rvorig, orient) == 0) return 0;
    forig = ((double) orient) * gv->dd(i, j)[family];

    if (consistency(rvorig, rvside, orient) == 0) return 0;
    fside = ((double) orient) * gv->dd(i + 1, j)[family];

    if (consistency(rvside, rvtop, orient) == 0) return 0;
    ftop = ((double) orient) * gv->dd(i, j + 1)[family];

    if (orient == -1) return 0;

    foncub[1] = forig; // Was: foncub[0][1]
    foncub[0] = fside; // Was: foncub[0][0]
    foncub[3] = ftop;  // Was: foncub[0][3]

    // Only follow for the fourth point if is square
    if (is_square == CELL_IS_TRIANGLE) return 1;

    // Test for real eigenvalues on the square
    if (!gv->eig_is_real(i + 1, j + 1)[family]) return 0;
    rvopps[0] = gv->e(i + 1, j + 1)[family].vrr[0];
    rvopps[1] = gv->e(i + 1, j + 1)[family].vrr[1];

    // Make vectors consistent at pairs of corners (if possible)
    // if vectors consistent, compute directional derivatives
    if (consistency(rvtop, rvopps, orient) == 0) return 0;
    fopps = ((double) orient) * gv->dd(i + 1, j + 1)[family];
    if (consistency(rvopps, rvside, orient) == 0) return 0;

    if (orient == -1) return 0;

    foncub[2] = fopps; // Was: foncub[0][2]

    return 1;
}

int Inflection_Curve::consistency(double *v1, double *v2, int &orient) {

    // When v2 is in opposite direction, must be changed. Orientation is negative
    if ((v1[0] * v2[0] + v1[1] * v2[1]) < 0.0) {
        v2[0] = -v2[0];
        v2[1] = -v2[1];
        orient = -1;
    } else {
        orient = 1;
    }

    // The angle between vectors must be lower than 45degrees ( cos45 = 0.7071... )
    // We are assuming normal vectors.
    // TODO: O angulo de 45graus eh arbitrario
//    if ((v1[0] * v2[0] + v1[1] * v2[1]) < 0.7071) return 0;

    return 1;
}

int Inflection_Curve::curve(const FluxFunction *f, const AccumulationFunction *a, 
                            GridValues &g, int fam, std::vector<RealVector> &inflection_curve) {
                            
    inflection_curve.clear();
    
    g.fill_dirdrv_on_grid(f, a);printf("Inflection_Curve::curve\n");
    gv = &g;

    // family MUST be a member of Inflection_Curve
    family = fam;

    int info = ContourMethod::contour2d(this, inflection_curve);
    
    cout<<"Tamanho da inflexao " <<inflection_curve.size()<<endl;

    return info;
}


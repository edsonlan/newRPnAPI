#include "Elliptic_Boundary.h"

Elliptic_Boundary::~Elliptic_Boundary() {
}

int Elliptic_Boundary::function_on_square(double *foncub, int i, int j) {
    int is_square = gv->cell_type(i, j);
    /* Here we write the discriminant for non-trivial accumulation as:
     *    DISCRIMINAT = (f11*g22 - f22*g11)^2 + 4*g11*g22*(f12-g12)*(f21-g21)
     * where fij stands for df_i/dx_j, as well for gij.
     */
    double f_aux[4];
    Matrix<double> JF, JG;

    JF.resize(2, 2); 
    JG.resize(2, 2);

    for (int k = 0; k < 2; k++) {
        for (int l = 0; l < 2; l++) {
            if ( (is_square == CELL_IS_TRIANGLE) && (k+l == 2) ) continue;
            for (int r = 0; r < 2; r++) {
                for (int s = 0; s < 2; s++) {
                    JF(r, s) = gv->JF_on_grid(i + k, j + l)(r, s);
                    JG(r, s) = gv->JG_on_grid(i + k, j + l)(r, s);
                }
            }
            double trace = JF(0, 0)*JG(1, 1) - JF(1, 1)*JG(0, 0);
            f_aux[2*l + k] = trace * trace
                           + 4 * JG(0, 0)*JG(1, 1) * ( JF(0, 1) - JG(0, 1) ) * ( JF(1, 0) - JG(1, 0) );
        }
    }

    foncub[1] = f_aux[0]; // Was: foncub[0][1]
    foncub[0] = f_aux[1]; // Was: foncub[0][0]
    foncub[3] = f_aux[2]; // Was: foncub[0][2]

    // Only useful if the cell is a square.
    //
    if (is_square == CELL_IS_SQUARE) foncub[2] = f_aux[3]; // Was: foncub[0][2]

    return 1;
}

int Elliptic_Boundary::curve(const FluxFunction *f, const AccumulationFunction *a, 
                          GridValues &g, std::vector<RealVector> &elliptic_boundary) {

    ff = f;
    aa = a;

    g.fill_Jacobians_on_grid(f, a);
    gv = &g;

    elliptic_boundary.clear();

    int info = ContourMethod::contour2d(this, elliptic_boundary);

    return info;
}

void Elliptic_Boundary::map(const RealVector &r, double &f, RealVector &map_Jacobian) {
    int n = r.size();
    map_Jacobian.resize(n);

    double F[n], JF[n][n], G[n], JG[n][n], HF[n][n][n], HG[n][n][n];

    RealVector p(r);

    ff->fill_with_jet(n, p.components(), 2, F, &JF[0][0], &HF[0][0][0]);
    aa->fill_with_jet(n, p.components(), 2, G, &JG[0][0], &HG[0][0][0]);

    // See the definition in Elliptic_Boundary::function_on_square; this is not exactly the trace!
    double trace = JF[0][0]*JG[1][1] - JF[1][1]*JG[0][0];
    double G_aux = JG[0][0]*JG[1][1];

    double aux_01 = JF[0][1] - JG[0][1];
    double aux_10 = JF[1][0] - JG[1][0];

    f = trace * trace + 4.0 * G_aux * aux_01 * aux_10;

    map_Jacobian.component(0) = 2.0 * trace * (HF[0][0][0]*JG[1][1] + JF[0][0]*HG[1][1][0] - HF[1][1][0]*JG[0][0] - JF[0][0]*HG[1][1][0])
                              + 4.0 * (HG[0][0][0]*JG[1][1] + JG[0][0]*HG[1][1][0]) * aux_01 * aux_10
                              + 4.0 * G_aux * ( (HF[0][1][0] - HG[0][1][0]) * aux_10 + aux_01 * (HF[1][0][0] - HG[1][0][0]) );

    map_Jacobian.component(1) = 2.0 * trace * (HF[0][0][1]*JG[1][1] + JF[0][0]*HG[1][1][1] - HF[1][1][1]*JG[0][0] - JF[0][0]*HG[1][1][1])
                              + 4.0 * (HG[0][0][1]*JG[1][1] + JG[0][0]*HG[1][1][1]) * aux_01 * aux_10
                              + 4.0 * G_aux * ( (HF[0][1][1] - HG[0][1][1]) * aux_10 + aux_01 * (HF[1][0][1] - HG[1][0][1]) );

   return;
}

bool Elliptic_Boundary::improvable(void) {
    return true;
}

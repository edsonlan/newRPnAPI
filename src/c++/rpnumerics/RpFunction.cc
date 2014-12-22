/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RpFunction.cc
 **/

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RpFunction.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */
RpFunction::~RpFunction(void){}

void RpFunction::fill_with_jet(int n, const double *in, int degree, double *F, double *J, double *H) const {
    RealVector r(n);
    for (int i = 0; i < n; i++) r.component(i) = in[i];

    // Will this work? There is a const somewhere in fluxParams.
    //FluxParams fp(r);
    //flux_object->fluxParams(FluxParams(r)); // flux_object->fluxParams(fp);

//    WaveState state_c(r);
    JetMatrix c_jet(n);

    jet(WaveState(r), c_jet, degree);

    // Fill F
    if (F != 0) for (int i = 0; i < n; i++) F[i] = c_jet.get(i);

    // Fill J
    if (J != 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                J[i * n + j] = c_jet.get(i, j);
            }
        }
    }

    // Fill H
    if (H != 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    H[(i * n + j) * n + k] = c_jet.get(i, j, k);
                }
            }
        }
    }

    return;
}


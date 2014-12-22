#include "Quad2HugoniotFunction.h"

Quad2HugoniotFunction::Quad2HugoniotFunction(const RealVector &U, const Quad2FluxFunction & fluxFunction) : HugoniotFunctionClass(fluxFunction) {

    int n = U.size();

    Uref.resize(n);
    for (int i = 0; i < n; i++) Uref.component(i) = U.component(i);

    // TODO: The flux object must be initialized somehow (be it created here or outside, etc.)
    UrefJetMatrix.resize(n);
    WaveState u(Uref); // TODO: Check this.
    fluxFunction.jet(u, UrefJetMatrix, 1);

    // Find the eigenpairs of the Jacobian of the jet at Uref
    double A[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = UrefJetMatrix.get(i, j);
        }
    }
    Eigen::eig(n, &A[0][0], ve_uref);

    double epsilon = 1.0e-6; // TODO: Tolerance must be created within a class which will be used by everyone. Say, "class Tolerance", or something like that.

    if (fabs(ve_uref[0].i) < epsilon) Uref_is_elliptic = false;
    else Uref_is_elliptic = true;
}

Quad2HugoniotFunction::~Quad2HugoniotFunction() {
    //    delete stone;
}

void Quad2HugoniotFunction::setReferenceVector(const RealVector & refVec) {
    Quad2FluxFunction & fluxFunction = (Quad2FluxFunction &) getFluxFunction();
    Uref = refVec;
    // TODO: The flux object must be initialized somehow (be it created here or outside, etc.)
    UrefJetMatrix.resize(refVec.size());
    int n = refVec.size();
    WaveState u(refVec); // TODO: Check this.
    fluxFunction.jet(u, UrefJetMatrix, 1);

    // Find the eigenpairs of the Jacobian of the jet at Uref
    double A[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = UrefJetMatrix.get(i, j);
        }
    }
    Eigen::eig(n, &A[0][0], ve_uref);

    double epsilon = 1.0e-6; // TODO: Tolerance must be created within a class which will be used by everyone. Say, "class Tolerance", or something like that.

    if (fabs(ve_uref[0].i) < epsilon) Uref_is_elliptic = false;
    else Uref_is_elliptic = true;

    HugoniotFunctionClass::setReferenceVector(refVec);


}

double Quad2HugoniotFunction::HugoniotFunction(const RealVector & U) {

    Quad2FluxFunction & fluxFunction = (Quad2FluxFunction &) getFluxFunction();

    double u = U(0);
    double uref = Uref(0);

    double v = U(1);
    double vref = Uref(1);

    double du = u - uref;
    double dv = v - vref;


    double epsilon = 1.0e-6; // TODO: Tolerance must be created within a
    //class which will be used by everyone. Say, "class Tolerance", or something
    //like that.

    // Find F at the point
    double dh1, dh2;
    double df1, df2;

    WaveState wu(U); // TODO: Is this correct?
    JetMatrix UJetMatrix(U.size());

    if (fabs(du) + fabs(dv) >= epsilon) {
        fluxFunction.jet(wu, UJetMatrix, 0);

        dh1 = du;
        dh2 = dv;
        df1 = UJetMatrix.get(0) - UrefJetMatrix.get(0);
        df2 = UJetMatrix.get(1) - UrefJetMatrix.get(1);
    } else {
        fluxFunction.jet(wu, UJetMatrix, 1);

        dh1 = du;
        dh2 = dv;
        df1 = UrefJetMatrix.get(0, 0) * du + UrefJetMatrix.get(0, 1) * dv;
        df2 = UrefJetMatrix.get(1, 0) * du + UrefJetMatrix.get(1, 1) * dv;
    }

    double hugoniot = dh2 * df1 - dh1*df2;

    if (fabs(hugoniot) <= epsilon * (fabs(dh2 * df1) + fabs(dh1 * df2)))
        hugoniot = 0.0;

    return hugoniot;

}


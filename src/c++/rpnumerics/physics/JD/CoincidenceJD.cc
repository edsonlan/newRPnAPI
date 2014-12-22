#include "CoincidenceJD.h"

CoincidenceJD::CoincidenceJD(const JDFluxFunction *f, const JDAccumulationFunction *a) : flux(f), accum(a) {
}

CoincidenceJD::~CoincidenceJD(){
}

void CoincidenceJD::lambdas(const RealVector &u, double &lambda_s, double &lambda_e, double &lambda_diff) const {
    lambda_s = this->lambda_s(u);
    lambda_e = this->lambda_e(u);

    lambda_diff = this->lambda_diff(u);

    return;
}

double CoincidenceJD::lambda_s(const RealVector &p) const {
    // lambda_s = df_ds.
    //
    JetMatrix f(2);

    flux->jet(p, f, 1);

    return f.get(0, 0);
}

double CoincidenceJD::lambda_e(const RealVector &p) const {
    JetMatrix f(2);
    flux->jet(p, f, 0);

    double phi = f.get(0);
    double alpha_dot = flux->alpha_dot(p);
    double u = p(0);

    return (phi + alpha_dot)/u;
}

double CoincidenceJD::lambda_diff(const RealVector &p) const {
    return lambda_s(p) - lambda_e(p);
}

bool CoincidenceJD::extension_basis(const RealVector &u, double &fe, double &se) const {
    return true;
}


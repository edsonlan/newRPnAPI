#include "ICDOWCoincidence.h"

ICDOWCoincidence::ICDOWCoincidence(ICDOWHydrodynamics *hy, Parameter *phi) : Coincidence() {
    hydro  = hy;
    phi_parameter = phi;
}

ICDOWCoincidence::~ICDOWCoincidence(){
}

void ICDOWCoincidence::lambdas(const RealVector &p, double &lambda_s, double &lambda_e, double &lambda_diff) const {
    double sw = p(0);
    double u  = p(2);

    double phi = phi_parameter->value();

    JetMatrix fw_jet;
    hydro->water_fractional_flow(sw, 1, fw_jet);
    double fw      = fw_jet.get(0);
    double dfw_dsw = fw_jet.get(0, 0);

    lambda_e = (u/phi)*fw/sw;
    lambda_s = (u/phi)*dfw_dsw;

    lambda_diff = lambda_s - lambda_e;

    return;
}

double ICDOWCoincidence::lambda_s(const RealVector &p) const {
    double ls, le, d;
    lambdas(p, ls, le, d);

    return ls;
}

double ICDOWCoincidence::lambda_e(const RealVector &p) const {
    double ls, le, d;
    lambdas(p, ls, le, d);

    return le;
}

double ICDOWCoincidence::lambda_diff(const RealVector &p) const {
    double ls, le, d;
    lambdas(p, ls, le, d);

    return d;
}


bool ICDOWCoincidence::extension_basis(const RealVector &u, double &fe, double &se) const {
    return true;
}


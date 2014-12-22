#include "JDEvap_Extension.h"

JDEvap_Extension::JDEvap_Extension(const JDFluxFunction *f, const Coincidence *c) : Extension(), flux_(f), coincidence_(c) {
    // We assume that the component which is constant along the extension is the second one.
    //
    index_of_constant = 1;

    // We assume that the component which is NOT constant along the extension is the first one.
    //
    index_of_non_constant = 0;
}

JDEvap_Extension::~JDEvap_Extension(){
}

int JDEvap_Extension::extension(const RealVector &p, RealVector &ext_p){
    ext_p = p;

    double uminus = p(0);

    double lambda_e_minus = coincidence_->lambda_e(p);

    double uplus  = lambda_e_minus + 2.0 - uminus;

    ext_p(0) = uplus;

    return EXTENSION_OK;
}


#include "SubPhysics.h"

SubPhysics::SubPhysics(){
    flux_ = 0;
    accumulation_ = 0;
    boundary_ = 0;
    viscosity_matrix_ = 0;
    
    rarefactioncurve_ = 0;
    compositecurve_ = 0;
    hugoniotcontinuation_ = 0;
    shockcurve_ = 0;
    wavecurvefactory_ = 0;

    explicitbifurcationcurves_ = 0;

    coincidence_ = 0;
    coincidence_contour_ = 0;

    inflection_curve_ = 0;

    bifurcationcurve_ = 0;

    doublecontact_ = 0;
}

SubPhysics::~SubPhysics(){
}




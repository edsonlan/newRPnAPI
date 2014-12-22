#include "ThreePhaseFlowViscosity.h"
#include "ThreePhaseFlowSubPhysics.h"

ThreePhaseFlowViscosity::ThreePhaseFlowViscosity(ThreePhaseFlowSubPhysics *t): subphysics_(t){
}

ThreePhaseFlowViscosity::~ThreePhaseFlowViscosity(){
}

double ThreePhaseFlowViscosity::water_viscosity(const RealVector &p){
    return subphysics_->muw()->value();
}

double ThreePhaseFlowViscosity::oil_viscosity(const RealVector &p){
    return subphysics_->muo()->value();
}

double ThreePhaseFlowViscosity::gas_viscosity(const RealVector &p){
    return subphysics_->mug()->value();
}


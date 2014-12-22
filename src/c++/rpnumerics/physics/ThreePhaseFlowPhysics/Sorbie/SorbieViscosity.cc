#include "SorbieViscosity.h"
#include "SorbieSubPhysics.h"

SorbieViscosity::SorbieViscosity(SorbieSubPhysics *t): ThreePhaseFlowViscosity((ThreePhaseFlowSubPhysics*)t){
}

SorbieViscosity::~SorbieViscosity(){
}

int SorbieViscosity::gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet){
    mug_jet.resize(2, 1);

    mug_jet.set(0, subphysics_->mug()->value());

    mug_jet.set(0, 0, 0.0);
    mug_jet.set(0, 1, 0.0);

    mug_jet.set(0, 0, 0, 0.0);
    mug_jet.set(0, 0, 1, 0.0);
    mug_jet.set(0, 1, 0, 0.0);
    mug_jet.set(0, 1, 1, 0.0);

    return VISCOSITY_OK;
}


#include "ThreePhaseFlowMobility.h"
#include "ThreePhaseFlowSubPhysics.h"

ThreePhaseFlowMobility::ThreePhaseFlowMobility(ThreePhaseFlowSubPhysics *s): permeability(s->permeability()), viscosity(s->viscosity()){
    subphysics = s;
}

ThreePhaseFlowMobility::~ThreePhaseFlowMobility(){
}

double ThreePhaseFlowMobility::water_mobility(const RealVector &p){
    JetMatrix perm;
    permeability->PermeabilityWater_jet(p, 0, perm);

    double visc = viscosity->water_viscosity(p);

    return perm.get(0)/visc;
}

double ThreePhaseFlowMobility::oil_mobility(const RealVector &p){
    JetMatrix perm;
    permeability->PermeabilityOil_jet(p, 0, perm);

    double visc = viscosity->oil_viscosity(p);

    return perm.get(0)/visc;
}

double ThreePhaseFlowMobility::gas_mobility(const RealVector &p){
    JetMatrix perm;
    permeability->PermeabilityGas_jet(p, 0, perm);

    double visc = viscosity->gas_viscosity(p);

    return perm.get(0)/visc;
}


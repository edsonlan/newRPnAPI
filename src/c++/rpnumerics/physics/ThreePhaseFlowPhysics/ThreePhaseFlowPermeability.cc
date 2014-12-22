#include "ThreePhaseFlowPermeability.h"
#include "ThreePhaseFlowSubPhysics.h"

ThreePhaseFlowPermeability::ThreePhaseFlowPermeability(ThreePhaseFlowSubPhysics *s): AuxiliaryFunction(),
        subphysics_(s){
}

ThreePhaseFlowPermeability::~ThreePhaseFlowPermeability(){
}


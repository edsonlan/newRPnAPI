#include "Flux2Comp2PhasesAdimensionalized_Params.h"

Flux2Comp2PhasesAdimensionalized_Params::Flux2Comp2PhasesAdimensionalized_Params(const RealVector & paramsVector, Thermodynamics* TD) : FluxParams(paramsVector),const_gravity(9.8) {
    TD_ = TD;
    if (paramsVector.component(2) == 1.0) {
        has_gravity_ = true;
    }
    else {
        has_gravity_ = false;
    }

    if (paramsVector.component(3) == 1.0) {
        has_horizontal_ = true;
    }
    else {
        has_horizontal_ = false;
    }


}

Flux2Comp2PhasesAdimensionalized_Params::~Flux2Comp2PhasesAdimensionalized_Params() {

}

Flux2Comp2PhasesAdimensionalized_Params::Flux2Comp2PhasesAdimensionalized_Params(const Flux2Comp2PhasesAdimensionalized_Params & copy) : FluxParams(copy.params()),const_gravity(9.8) {
    component(0, copy.abs_perm);
    component(1, copy.sin_beta);
    TD_ = copy.TD_;
    has_gravity_ = copy.has_gravity_;
    has_horizontal_ = copy.has_horizontal_;
}

Thermodynamics * Flux2Comp2PhasesAdimensionalized_Params::get_thermodynamics(void) const {
    return TD_;
}

bool Flux2Comp2PhasesAdimensionalized_Params::has_gravity(void) const {
    return has_gravity_;
}

bool Flux2Comp2PhasesAdimensionalized_Params::has_horizontal(void) const {
    return has_horizontal_;
}


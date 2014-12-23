#include "SorbiePermeability.h"

SorbiePermeability::SorbiePermeability(ThreePhaseFlowSubPhysics *s): ThreePhaseFlowPermeability(s){
}

SorbiePermeability::~SorbiePermeability(){
}

// The reduced permeabilities must be changed after changing this method!
//
int SorbiePermeability::PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water){
    water.resize(2, 1);

    double sw = state(0);
    double so = state(1);
    double sg = 1.0 - sw - so;

        if (degree >= 0){
            double kw = sw*sg;
            water.set(0, kw);

            if (degree >= 1){
                double dkw_dsw = 1.0 - 2.0*sw - so;
                double dkw_dso = -sw;

                water.set(0, 0, dkw_dsw);
                water.set(0, 1, dkw_dso);

                if (degree >= 2){
                    double d2kw_dsw2  = -2.0;
                    double d2kw_dswso = -1.0;
                    double d2kw_dsosw = -1.0;
                    double d2kw_dso2  =  0.0;

                    water.set(0, 0, 0, d2kw_dsw2);
                    water.set(0, 0, 1, d2kw_dswso);
                    water.set(0, 1, 0, d2kw_dsosw);
                    water.set(0, 1, 1, d2kw_dso2);
                }
            }
        }

    return degree;
}

// The reduced permeabilities must be changed after changing this method!
//
int SorbiePermeability::PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil){
    oil.resize(2, 1);

//    double sw = state(0);
    double so = state(1);
//    double sg = 1.0 - sw - so;

        if (degree >= 0){
            double ko = so*so;
            oil.set(0, ko);

            if (degree >= 1){
                oil.set(0, 0, 0.0);    // dkw_dsw
                oil.set(0, 1, 2.0*so); // dkw_dso

                if (degree >= 2){
                    oil.set(0, 0, 0, 0.0); // d2kw_dsw2
                    oil.set(0, 0, 1, 0.0); // d2kw_dswso
                    oil.set(0, 1, 0, 0.0); // d2kw_dsosw
                    oil.set(0, 1, 1, 2.0); // d2kw_dso2
                }
            }
        }

    return degree;
}

// The reduced permeabilities must be changed after changing this method!
//
int SorbiePermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas){
    gas.resize(2, 1);

    double sw = state(0);
    double so = state(1);
    double sg = 1.0 - sw - so;

        if (degree >= 0){
            double kg = sg*sg;
            gas.set(0, kg);

            if (degree >= 1){
                gas.set(0, 0, -2.0*sg); // dkg_dsw
                gas.set(0, 1, -2.0*sg); // dkg_dso

                if (degree >= 2){
                    gas.set(0, 0, 0, 2.0); // d2kg_dsw2
                    gas.set(0, 0, 1, 2.0); // d2kg_dswso
                    gas.set(0, 1, 0, 2.0); // d2kg_dsosw
                    gas.set(0, 1, 1, 2.0); // d2kg_dso2
                }
            }
        }

    return degree;
}

void SorbiePermeability::reduced_permeability(const RealVector &state, RealVector &reduced){
    double rw, ro, rg;

    JetMatrix water, oil, gas;
    PermeabilityWater_jet(state, 0, water);
    PermeabilityGas_jet(state, 0, gas);
    PermeabilityOil_jet(state, 0, oil);

//    reduced(0) = water.get(0)/state(0);
//    reduced(1) = oil.get(0)/state(1);
//    reduced(2) = gas.get(0)/(1.0 - state(0) - state(1));

    double sw = state(0);
    double so = state(1);
    double sg = 1.0 - state(0) - state(1);

    reduced.resize(3);
    reduced(0) = sg;
    reduced(1) = so;
    reduced(2) = sg;

    return;
}


#ifndef _COREYQUADPERMEABILITY_
#define _COREYQUADPERMEABILITY_

#include "ThreePhaseFlowPermeability.h"

class CoreyQuadSubPhysics;

class CoreyQuadPermeability: public ThreePhaseFlowPermeability {
    private:
    protected:
    public:
        CoreyQuadPermeability(ThreePhaseFlowSubPhysics *s);
        virtual ~CoreyQuadPermeability();

        inline int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water);
        inline int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &water);
        inline int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &water);

        inline void reduced_permeability(const RealVector &state, RealVector &rp);
};

int CoreyQuadPermeability::PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water){
    water.resize(2, 1);

    double sw = state(0);
//    double so = state(1);
//    double sg = 1.0 - sw - so;

        if (degree >= 0){
            double kw = sw*sw;
            water.set(0, kw);

            if (degree >= 1){
                double dkw_dsw = 2.0*sw;
                double dkw_dso = 0.;

                water.set(0, 0, dkw_dsw);
                water.set(0, 1, dkw_dso);

                if (degree >= 2){
                    double d2kw_dsw2  = 2.0;
                    double d2kw_dswso = 0.0;
                    double d2kw_dsosw = 0.0;
                    double d2kw_dso2  = 0.0;

                    water.set(0, 0, 0, d2kw_dsw2);
                    water.set(0, 0, 1, d2kw_dswso);
                    water.set(0, 1, 0, d2kw_dsosw);
                    water.set(0, 1, 1, d2kw_dso2);
                }
            }
        }

    return degree;
}

int CoreyQuadPermeability::PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil){
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

int CoreyQuadPermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas){
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

void CoreyQuadPermeability::reduced_permeability(const RealVector &state, RealVector &rp){
    double sw = state(0);
    double so = state(1);
    double sg = 1.0 - sw - so;

    rp.resize(3);

    // Water.
    //
    rp(0) = sw;

    // Oil.
    //
    rp(1) = so;

    // Gas.
    //
    rp(2) = sg;

    return;
}

#endif // _COREYQUADPERMEABILITY_


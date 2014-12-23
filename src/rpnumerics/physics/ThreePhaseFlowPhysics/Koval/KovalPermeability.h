#ifndef _KOVALPERMEABILITY_
#define _KOVALPERMEABILITY_

#include "ThreePhaseFlowPermeability.h"

class KovalSubPhysics;

class KovalPermeability: public ThreePhaseFlowPermeability {
    private:
    protected:
    public:
        KovalPermeability(ThreePhaseFlowSubPhysics *s);
        virtual ~KovalPermeability();

        inline int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water);
        inline int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &water);
        inline int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &water);

        inline void reduced_permeability(const RealVector &state, RealVector &rp);
};

int KovalPermeability::PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water){
    water.resize(1);
    double sw = state(0);
//    double so = state(1);
//    double sg = 1.0 - sw - so;

        if (degree >= 0){
            double kw = sw;
            water.set(0, kw);

            if (degree >= 1){
                double dkw_dsw = 1.0;
                double dkw_dso = 0.;

                water.set(0, 0, dkw_dsw);
                water.set(0, 1, dkw_dso);

                if (degree >= 2){
                    double d2kw_dsw2  = 0.0;
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

int KovalPermeability::PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil){
    oil.resize(1);
//    double sw = state(0);
    double so = state(1);
//    double sg = 1.0 - sw - so;

        if (degree >= 0){
            double ko = so;
            oil.set(0, ko);

            if (degree >= 1){
                oil.set(0, 0, 0.0); // dko_dsw
                oil.set(0, 1, 1.0); // dko_dso

                if (degree >= 2){
                    oil.set(0, 0, 0, 0.0); // d2ko_dsw2
                    oil.set(0, 0, 1, 0.0); // d2ko_dswso
                    oil.set(0, 1, 0, 0.0); // d2ko_dsosw
                    oil.set(0, 1, 1, 0.0); // d2ko_dso2
                }
            }
        }

    return degree;
}

int KovalPermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas){
    gas.resize(1);

    double sw = state(0);
    double so = state(1);
    double sg = 1.0 - sw - so;

        if (degree >= 0){
            double kg = sg;
            gas.set(0, kg);

            if (degree >= 1){
                gas.set(0, 0, -1.0); // dkg_dsw
                gas.set(0, 1, -1.0); // dkg_dso

                if (degree >= 2){
                    gas.set(0, 0, 0, 0.0); // d2kg_dsw2
                    gas.set(0, 0, 1, 0.0); // d2kg_dswso
                    gas.set(0, 1, 0, 0.0); // d2kg_dsosw
                    gas.set(0, 1, 1, 0.0); // d2kg_dso2
                }
            }
        }

    return degree;
}

void KovalPermeability::reduced_permeability(const RealVector &state, RealVector &rp){
    rp.resize(3);

    // Water.
    //
    rp(0) = 1.0;

    // Oil.
    //
    rp(1) = 1.0;

    // Gas.
    //
    rp(2) = 1.0;

    return;
}

#endif // _KOVALPERMEABILITY_


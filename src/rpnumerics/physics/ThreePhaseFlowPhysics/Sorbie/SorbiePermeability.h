#ifndef _SORBIEPERMEABILITY_
#define _SORBIEPERMEABILITY_

#include "ThreePhaseFlowPermeability.h"

class SorbieSubPhysics;

class SorbiePermeability: public ThreePhaseFlowPermeability {
    private:
    protected:
    public:
        SorbiePermeability(ThreePhaseFlowSubPhysics *s);
        virtual ~SorbiePermeability();

        inline int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water);
        inline int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &water);
        inline int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &water);

        inline void reduced_permeability(const RealVector &state, RealVector &rp);
};

#endif // _SORBIEPERMEABILITY_


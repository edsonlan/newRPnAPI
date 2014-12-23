#ifndef _STONEVISCOSITY_
#define _STONEVISCOSITY_

#include "ThreePhaseFlowViscosity.h"

class StoneSubPhysics;

class StoneViscosity: public ThreePhaseFlowViscosity {
    private:
    protected:
    public:
        StoneViscosity(StoneSubPhysics *t);
        virtual ~StoneViscosity();

        int gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet);
};

#endif // _STONEVISCOSITY_


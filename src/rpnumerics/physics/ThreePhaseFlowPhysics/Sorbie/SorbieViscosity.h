#ifndef _SORBIEVISCOSITY_
#define _SORBIEVISCOSITY_

#include "ThreePhaseFlowViscosity.h"

class SorbieSubPhysics;

class SorbieViscosity: public ThreePhaseFlowViscosity {
    private:
    protected:
    public:
        SorbieViscosity(SorbieSubPhysics *t);
        virtual ~SorbieViscosity();

        int gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet);
};

#endif // _SORBIEVISCOSITY_


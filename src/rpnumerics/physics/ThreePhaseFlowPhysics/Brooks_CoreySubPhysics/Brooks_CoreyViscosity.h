#ifndef _BROOKS_COREYVISCOSITY_
#define _BROOKS_COREYVISCOSITY_

#include "ThreePhaseFlowViscosity.h"

class Brooks_CoreySubPhysics;

class Brooks_CoreyViscosity: public ThreePhaseFlowViscosity {
    private:
    protected:
    public:
        Brooks_CoreyViscosity(Brooks_CoreySubPhysics *t);
        virtual ~Brooks_CoreyViscosity();

        int gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet);
};

#endif // _BROOKS_COREYVISCOSITY_


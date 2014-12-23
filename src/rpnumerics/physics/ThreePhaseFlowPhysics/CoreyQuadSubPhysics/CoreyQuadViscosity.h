#ifndef _COREYQUADVISCOSITY_
#define _COREYQUADVISCOSITY_

#include "ThreePhaseFlowViscosity.h"

class CoreyQuadSubPhysics;

class CoreyQuadViscosity: public ThreePhaseFlowViscosity {
    private:
    protected:
    public:
        CoreyQuadViscosity(CoreyQuadSubPhysics *t);
        virtual ~CoreyQuadViscosity();

        int gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet);
};

#endif // _COREYQUADVISCOSITY_


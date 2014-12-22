#ifndef _KOVALVISCOSITY_
#define _KOVALVISCOSITY_

#include "ThreePhaseFlowViscosity.h"

class KovalSubPhysics;

class KovalViscosity: public ThreePhaseFlowViscosity {
    private:
    protected:
    public:
        KovalViscosity(KovalSubPhysics *t);
        virtual ~KovalViscosity();

        int gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet);
};

#endif // _KOVALVISCOSITY_


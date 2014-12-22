#ifndef _STONESUBPHYSICS_
#define _STONESUBPHYSICS_

#include "ThreePhaseFlowSubPhysics.h"
#include "StoneFluxFunction.h"
#include "Stone_Explicit_Bifurcation_Curves.h"
#include "LSODE.h"
#include "StoneViscosity.h"

class StoneSubPhysics : public ThreePhaseFlowSubPhysics {
    private:
    protected:
//        StonePermeability *permeability_;
    public:
        StoneSubPhysics();
        ~StoneSubPhysics();
};

#endif // _STONESUBPHYSICS_


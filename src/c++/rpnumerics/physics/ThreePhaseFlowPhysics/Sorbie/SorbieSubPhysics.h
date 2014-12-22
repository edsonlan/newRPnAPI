#ifndef _SORBIESUBPHYSICS_
#define _SORBIESUBPHYSICS_

#include "ThreePhaseFlowSubPhysics.h"
#include "SorbieFluxFunction.h"
#include "LSODE.h"
#include "SorbiePermeability.h"
#include "SorbieViscosity.h"
#include "WaveCurveFactory.h" 

class SorbieSubPhysics : public ThreePhaseFlowSubPhysics {
    private:
    protected:
    public:
        SorbieSubPhysics();
        virtual ~SorbieSubPhysics();
};

#endif // _SORBIESUBPHYSICS_


#ifndef _KOVALSUBPHYSICS_
#define _KOVALSUBPHYSICS_

#include "ThreePhaseFlowSubPhysics.h"
#include "LSODE.h"
#include "KovalPermeability.h"

#include "KovalFluxFunction.h"
#include "KovalViscosity.h"

class KovalSubPhysics : public ThreePhaseFlowSubPhysics {
    private:
    protected:
    public:
        KovalSubPhysics();
        virtual ~KovalSubPhysics();
};

#endif // _KOVALSUBPHYSICS_


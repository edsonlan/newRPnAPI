#ifndef _BROOKS_COREYSUBPHYSICS_
#define _BROOKS_COREYSUBPHYSICS_

#include "ThreePhaseFlowSubPhysics.h"
#include "Brooks_CoreyPermeability.h"
#include "Brooks_CoreyFluxFunction.h"
#include "LSODE.h"
#include "Brooks_CoreyViscosity.h"

#define BROOKS_COREYGENERICPOINT 0

class Brooks_CoreySubPhysics : public ThreePhaseFlowSubPhysics {
    private:
    protected:
        Parameter *lambda_parameter;
        Parameter *cnw_parameter, *cno_parameter, *cng_parameter;

    public:
        Brooks_CoreySubPhysics();
        virtual ~Brooks_CoreySubPhysics();

        void shock_cases(std::vector<int> &type, std::vector<std::string> &name) const;
};

#endif // _BROOKS_COREYSUBPHYSICS_


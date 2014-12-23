#ifndef _THREEPHASEFLOWFLUXFUNCTION_
#define _THREEPHASEFLOWFLUXFUNCTION_

#include "Parameter.h"
#include "FluxFunction.h"

class ThreePhaseFlowFluxFunction : public FluxFunction {
    private:
    protected:
        Parameter *muw_, *muo_, *mug_;

    public:
        ThreePhaseFlowFluxFunction();
        virtual ~ThreePhaseFlowFluxFunction();

        virtual Parameter *muw(){return muw_;}
        virtual Parameter *muo(){return muo_;}
        virtual Parameter *mug(){return mug_;}
};

#endif // _THREEPHASEFLOWFLUXFUNCTION_


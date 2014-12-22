#ifndef _BROOKS_COREYFLUXFUNCTION_
#define _BROOKS_COREYFLUXFUNCTION_

#include "FluxFunction.h"
#include "Brooks_CoreyPermeability.h"

class Brooks_CoreyFluxFunction : public FluxFunction {
    private:
    protected:
        Parameter *muw_parameter, *muo_parameter, *mug_parameter;
        Brooks_CoreyPermeability *perm;
    public:
        Brooks_CoreyFluxFunction(Parameter *muw, Parameter *muo, Parameter *mug, Brooks_CoreyPermeability *p);
        virtual ~Brooks_CoreyFluxFunction();

        int jet(const WaveState &u, JetMatrix &f, int degree) const;

};

#endif // _BROOKS_COREYFLUXFUNCTION_


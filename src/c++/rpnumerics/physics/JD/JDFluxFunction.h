#ifndef _JDFLUXFUNCTION_
#define _JDFLUXFUNCTION_

#include "FluxFunction.h"
#include "Parameter.h"
#include "ImplicitHugoniotCurve.h"

class JDFluxFunction : public FluxFunction {
    private:
    protected:
        Parameter *epsilon_parameter;
    public:
        JDFluxFunction(Parameter *e);

        virtual ~JDFluxFunction();

        int jet(const WaveState &w, JetMatrix &f, int degree) const;
        JDFluxFunction * clone() const;

        // For the coincidence:
        //
        double alpha_dot(const RealVector &p) const;
};

#endif // _JDFLUXFUNCTION_


#ifndef _QUAD2FLUXFUNCTION_
#define _QUAD2FLUXFUNCTION_

#include "FluxFunction.h"
#include "Parameter.h"

class Quad2FluxFunction: public FluxFunction {
    private:
    protected:
        Parameter *a1_parameter, *b1_parameter, *c1_parameter, *d1_parameter, *e1_parameter;
        Parameter *a2_parameter, *b2_parameter, *c2_parameter, *d2_parameter, *e2_parameter;
    public:
        Quad2FluxFunction(Parameter *a1, Parameter *b1, Parameter *c1, Parameter *d1, Parameter *e1,
                          Parameter *a2, Parameter *b2, Parameter *c2, Parameter *d2, Parameter *e2);
        virtual ~Quad2FluxFunction();

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _QUAD2FLUXFUNCTION_


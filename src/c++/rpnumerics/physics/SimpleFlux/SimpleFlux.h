#ifndef _SIMPLEFLUX_
#define _SIMPLEFLUX_

#include "FluxFunction.h"

class SimpleFlux : public FluxFunction{
    protected:
        RealVector coef, coef_der1, coef_der2;
        void init(int n, const double *c);
    public:
        SimpleFlux(const RealVector &c);
        SimpleFlux(int n, const double *c);

        ~SimpleFlux();

        SimpleFlux * clone() const;

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _SIMPLEFLUX_


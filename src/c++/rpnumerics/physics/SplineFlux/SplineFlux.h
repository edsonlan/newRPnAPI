#ifndef _SPLINEFLUX_
#define _SPLINEFLUX_

#include "spline1d.h"
#include "FluxFunction.h"

class SplineFlux : public FluxFunction {
    private:
    protected:
        spline1dinterpolant spline;

        std::vector<RealVector> points;
    public:
        SplineFlux(const std::vector<RealVector> &p);
        virtual ~SplineFlux();

        SplineFlux * clone() const;

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _SPLINEFLUX_


#ifndef _IONFLUX_
#define _IONFLUX_

#include "FluxFunction.h"

#include "IonRatios.h"
#include "IonPermeability.h"

class IonFlux : public FluxFunction {
    private:
    protected:
        const IonRatios       *ratios;
        const IonPermeability *permeability;
    public:
        IonFlux(const IonRatios *r, const IonPermeability *p);
        virtual ~IonFlux();

        int jet(const WaveState &w, JetMatrix &m, int degree) const;

        IonFlux * clone() const;
};

#endif // _IONFLUX_


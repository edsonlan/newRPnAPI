#ifndef _DEADVOLATILEVOLATILEGASFLUXFUNCTION_
#define _DEADVOLATILEVOLATILEGASFLUXFUNCTION_

#include "FluxFunction.h"
#include "DeadVolatileVolatileGasThermodynamics.h"
#include "DeadVolatileVolatileGasHydrodynamics.h"

class DeadVolatileVolatileGasFluxFunction : public FluxFunction {
    private:
    protected:
        DeadVolatileVolatileGasThermodynamics *thermo;
        DeadVolatileVolatileGasHydrodynamics  *hydro;
    public:
        DeadVolatileVolatileGasFluxFunction(DeadVolatileVolatileGasThermodynamics *th, DeadVolatileVolatileGasHydrodynamics *hy);
        virtual ~DeadVolatileVolatileGasFluxFunction();

        int reduced_jet(const WaveState &u, JetMatrix &m, int degree) const;
        int reduced_jet(const RealVector &u, JetMatrix &m, int degree) const {
            WaveState w(u);

            int info = reduced_jet(w, m, degree);

            return info;
        }

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _DEADVOLATILEVOLATILEGASFLUXFUNCTION_


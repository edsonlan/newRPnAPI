#ifndef _DEADVOLATILEVOLATILEGASHYDRODYNAMICS_
#define _DEADVOLATILEVOLATILEGASHYDRODYNAMICS_

#include "AuxiliaryFunction.h"
#include "DeadVolatileVolatileGasThermodynamics.h"

class DeadVolatileVolatileGasHydrodynamics : public AuxiliaryFunction {
    private:
    protected:
        DeadVolatileVolatileGasThermodynamics *thermo;
    public:
        DeadVolatileVolatileGasHydrodynamics(DeadVolatileVolatileGasThermodynamics *t);
        virtual ~DeadVolatileVolatileGasHydrodynamics();

        void oil_fractional_flow(int degree, double s, double r, JetMatrix &fo);
        void fractional_flow(int degree, double s, double y, JetMatrix &f_jet);
};

#endif // _DEADVOLATILEVOLATILEGASHYDRODYNAMICS_


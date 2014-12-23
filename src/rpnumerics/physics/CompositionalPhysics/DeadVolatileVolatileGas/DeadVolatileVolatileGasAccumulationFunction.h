#ifndef _DEADVOLATILEVOLATILEGASACCUMULATIONFUNCTION_
#define _DEADVOLATILEVOLATILEGASACCUMULATIONFUNCTION_

#include "AccumulationFunction.h"
#include "DeadVolatileVolatileGasThermodynamics.h"

class DeadVolatileVolatileGasAccumulationFunction : public AccumulationFunction {
    private:
    protected:
        DeadVolatileVolatileGasThermodynamics *thermo;
        
        Parameter *phi_parameter;
    public:
        DeadVolatileVolatileGasAccumulationFunction(Parameter *phi, DeadVolatileVolatileGasThermodynamics *th);
        virtual ~DeadVolatileVolatileGasAccumulationFunction();

        int reduced_jet(const WaveState &state, JetMatrix &m, int degree) const;
        int reduced_jet(const RealVector &u, JetMatrix &m, int degree) const {
            WaveState w(u);

            int info = reduced_jet(w, m, degree);

            return info;
        }

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _DEADVOLATILEVOLATILEGASACCUMULATIONFUNCTION_


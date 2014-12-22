#ifndef _ACCUMULATIONSINGLEPHASELIQUIDADIMENSIONALIZED_
#define _ACCUMULATIONSINGLEPHASELIQUIDADIMENSIONALIZED_

#include <stdio.h>
#include <stdlib.h>
#include "AccumulationFunction.h"


#include "Thermodynamics.h"

class AccumulationSinglePhaseLiquidAdimensionalized : public AccumulationFunction {
    private:


        // Thermodynamics
        Thermodynamics *thermo;
    protected:
        double phi_;
    public:
//        AccumulationSinglePhaseLiquidAdimensionalized(const AccumulationSinglePhaseLiquidAdimensionalized &);
        AccumulationSinglePhaseLiquidAdimensionalized(double phi, Thermodynamics *t);
        AccumulationSinglePhaseLiquidAdimensionalized * clone() const;

        ~AccumulationSinglePhaseLiquidAdimensionalized();

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _ACCUMULATIONSINGLEPHASELIQUIDADIMENSIONALIZED_


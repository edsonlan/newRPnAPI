#ifndef _ACCUMULATIONSINGLEPHASEVAPORADIMENSIONALIZED_
#define _ACCUMULATIONSINGLEPHASEVAPORADIMENSIONALIZED_

#include <stdio.h>
#include <stdlib.h>
#include "AccumulationFunction.h"


#include "Thermodynamics.h"

class AccumulationSinglePhaseVaporAdimensionalized : public AccumulationFunction {
    private:
        // Thermodynamics
        Thermodynamics *thermo;
    protected:
        double phi_;
    public:
//        AccumulationSinglePhaseVaporAdimensionalized(const AccumulationSinglePhaseVaporAdimensionalized &);
        AccumulationSinglePhaseVaporAdimensionalized(double phi, Thermodynamics *t);
        AccumulationSinglePhaseVaporAdimensionalized * clone() const;

        ~AccumulationSinglePhaseVaporAdimensionalized();

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _ACCUMULATIONSINGLEPHASEVAPORADIMENSIONALIZED_


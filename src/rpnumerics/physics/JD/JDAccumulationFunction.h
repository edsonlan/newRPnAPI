#ifndef _JDACCUMULATIONFUNCTION_
#define _JDACCUMULATIONFUNCTION_

#include "AccumulationFunction.h"

class JDAccumulationFunction : public AccumulationFunction {
    private:
    protected:
    public:
        JDAccumulationFunction();

        virtual ~JDAccumulationFunction();

        int jet(const WaveState &w, JetMatrix &f, int degree) const;
        JDAccumulationFunction * clone() const;
};

#endif // _JDACCUMULATIONFUNCTION_


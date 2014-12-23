#ifndef _TRIVIALACCUMULATIONFUNCTION_
#define _TRIVIALACCUMULATIONFUNCTION_

#include "AccumulationFunction.h"

// This class provides a trivial accumulation. It is 
// expected to be a final class.
//
class TrivialAccumulationFunction : public AccumulationFunction {
    private:
    protected:
    public:
        TrivialAccumulationFunction();
        ~TrivialAccumulationFunction();

        int jet(const WaveState &w, JetMatrix &a, int degree) const;
        TrivialAccumulationFunction *clone() const;
};

#endif // _TRIVIALACCUMULATIONFUNCTION_


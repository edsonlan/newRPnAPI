#ifndef _QUAD2ACCUMULATIONFUNCTION_
#define	_QUAD2ACCUMULATIONFUNCTION_

#include "AccumulationFunction.h"

class Quad2AccumulationFunction: public AccumulationFunction {
    public:
        Quad2AccumulationFunction();
        virtual ~Quad2AccumulationFunction();
        
        int jet(const WaveState &w, JetMatrix &m, int degree) const;
};

#endif // _QUAD2ACCUMULATIONFUNCTION_


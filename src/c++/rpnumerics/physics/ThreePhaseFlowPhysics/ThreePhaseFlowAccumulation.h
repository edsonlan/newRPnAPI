#ifndef _THREEPHASEFLOWACCUMULATION_
#define _THREEPHASEFLOWACCUMULATION_

#include "AccumulationFunction.h"

class ThreePhaseFlowAccumulation : public AccumulationFunction {
    private:
    protected:
    public:
        ThreePhaseFlowAccumulation();
        virtual ~ThreePhaseFlowAccumulation();

        int jet(const WaveState &w, JetMatrix &m, int degree) const;
};

#endif // _THREEPHASEFLOWACCUMULATION_


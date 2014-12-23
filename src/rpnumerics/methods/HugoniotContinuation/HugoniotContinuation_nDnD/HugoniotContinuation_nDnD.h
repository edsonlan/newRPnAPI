#ifndef _HUGONIOTCONTINUATION_NDND_
#define _HUGONIOTCONTINUATION_NDND_

#include <iostream>

#include "HugoniotContinuation.h"

class HugoniotContinuation_nDnD : public HugoniotContinuation {
    private:
    protected:
        // Undefined in the base class:
        //
        //int fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &hyperplane, RealVector &Hugoniot_direction);

        RealVector orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane);
    public:
        HugoniotContinuation_nDnD(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb);
        virtual ~HugoniotContinuation_nDnD();

        // Undefined in the base class:
        //
        // void jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH); // Is no more.
        void jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                          const RealVector &G, const DoubleMatrix &JG,
                          const RealVector &C, const DoubleMatrix &JC,
                          RealVector &H, DoubleMatrix &nablaH);
};

#endif // _HUGONIOTCONTINUATION_NDND_

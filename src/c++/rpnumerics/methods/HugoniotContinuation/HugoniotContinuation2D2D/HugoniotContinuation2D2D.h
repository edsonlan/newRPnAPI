#ifndef _HUGONIOTCONTINUATION2D2D_
#define _HUGONIOTCONTINUATION2D2D_

#include <iostream>

#include "HugoniotContinuation_nDnD.h"

class HugoniotContinuation2D2D : public HugoniotContinuation_nDnD {
    private:
    protected:
        // Undefined in the base class:
        //
        //int fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &subspace, RealVector &Hugoniot_direction);
        RealVector orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane);
    public:
        HugoniotContinuation2D2D(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb);
        virtual ~HugoniotContinuation2D2D();

        // Undefined in the base class:
        //
        // void jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH); // Is no more.
        void jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                          const RealVector &G, const DoubleMatrix &JG,
                          const RealVector &C, const DoubleMatrix &JC,
                          RealVector &H, DoubleMatrix &nablaH);
};

#endif // _HUGONIOTCONTINUATION2D2D_


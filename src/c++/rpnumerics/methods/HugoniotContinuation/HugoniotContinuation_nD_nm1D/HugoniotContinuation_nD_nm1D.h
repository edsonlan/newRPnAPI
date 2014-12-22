#ifndef _HUGONIOTCONTINUATION_ND_NM1D_
#define HUGONIOTCONTINUATION_ND_NM1D

#include "HugoniotContinuation.h"

class HugoniotContinuation_nD_nm1D : public HugoniotContinuation {
    private:
    protected:
        // Undefined in the base class:
        //
        RealVector orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane);
    public:
        HugoniotContinuation_nD_nm1D(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb);
        virtual ~HugoniotContinuation_nD_nm1D();

        // Undefined in the base class:
        //
        void jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                          const RealVector &G, const DoubleMatrix &JG,
                          const RealVector &C, const DoubleMatrix &JC,
                          RealVector &H, DoubleMatrix &nablaH);
};

#endif // HUGONIOTCONTINUATION_ND_NM1D


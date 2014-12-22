#ifndef _HUGONIOTCONTINUATION3D2D_
#define _HUGONIOTCONTINUATION3D2D_

#include <iostream>

#include "HugoniotContinuation.h"

class HugoniotContinuation3D2D : public HugoniotContinuation {
    private:
    protected:
        // Undefined in the base class:
        //
        //int fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &hyperplane, RealVector &Hugoniot_direction);
        RealVector orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane);
    public:
        HugoniotContinuation3D2D(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb);
        virtual ~HugoniotContinuation3D2D();

        // Undefined in the base class:
        //
        // void jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH); // Is no more.
        void jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                          const RealVector &G, const DoubleMatrix &JG,
                          const RealVector &C, const DoubleMatrix &JC,
                          RealVector &H, DoubleMatrix &nablaH);

};

#endif // _HUGONIOTCONTINUATION3D2D_

#ifndef _THREEPHASEFLOWHUGONIOTCONTINUATION_
#define _THREEPHASEFLOWHUGONIOTCONTINUATION_

#include "HugoniotContinuation.h"

class ThreePhaseFlowSubPhysics;

class ThreePhaseFlowHugoniotContinuation: public HugoniotContinuation {
    private:
    protected:
        ThreePhaseFlowSubPhysics *subphysics;

        void jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                          const RealVector &G, const DoubleMatrix &JG, 
                          const RealVector &C, const DoubleMatrix &JC, 
                          RealVector &H, DoubleMatrix &nablaH);

        virtual RealVector orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane);

        static double H_water_vertex(ThreePhaseFlowHugoniotContinuation *obj, const RealVector &p);
        static double H_oil_vertex(ThreePhaseFlowHugoniotContinuation *obj, const RealVector &p);
        static double H_gas_vertex(ThreePhaseFlowHugoniotContinuation *obj, const RealVector &p);
    public:
        ThreePhaseFlowHugoniotContinuation(ThreePhaseFlowSubPhysics *s);
        virtual ~ThreePhaseFlowHugoniotContinuation();
};

#endif // _THREEPHASEFLOWHUGONIOTCONTINUATION_


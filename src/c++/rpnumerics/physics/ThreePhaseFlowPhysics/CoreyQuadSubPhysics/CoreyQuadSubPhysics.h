#ifndef _COREYQUADSUBPHYSICS_
#define _COREYQUADSUBPHYSICS_

#include "ThreePhaseFlowSubPhysics.h"
#include "CoreyQuadFluxFunction.h"
#include "LSODE.h"
#include "CoreyQuadExplicitHugoniotCurve.h"
#include "CoreyQuadPermeability.h"
#include "CoreyQuadTransitionalLine.h"
#include "CoreyQuadViscosity.h"
#include "CoreyQuadWaveCurveFactory.h" 

class CoreyQuadSubPhysics : public ThreePhaseFlowSubPhysics {
    private:
    protected:
    public:
        CoreyQuadSubPhysics();
        virtual ~CoreyQuadSubPhysics();

        // These three methods below will, in all probability, move to CoreyQuadSubPhysics.
        //
        RealVector E(){
            RealVector Ep(2);
            Ep(0) = 0.0;
            Ep(1) = muo_parameter->value()/(muo_parameter->value() + mug_parameter->value());

            return Ep;
        }

        RealVector B(){
            RealVector Bp(2);
            Bp(0) = muw_parameter->value()/(muw_parameter->value() + mug_parameter->value());
            Bp(1) = 0.0;

            return Bp;
        }

        RealVector D(){
            RealVector Dp(2);
            Dp(0) = muw_parameter->value()/(muw_parameter->value() + muo_parameter->value());
            Dp(1) = muo_parameter->value()/(muw_parameter->value() + muo_parameter->value());

            return Dp;
        }
};

#endif // _COREYQUADSUBPHYSICS_


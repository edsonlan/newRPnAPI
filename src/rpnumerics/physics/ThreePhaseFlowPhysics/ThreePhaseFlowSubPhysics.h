#ifndef _THREEPHASEFLOWSUBPHYSICS_
#define _THREEPHASEFLOWSUBPHYSICS_

#include "SubPhysics.h"

#include "IsoTriang2DBoundary.h"
#include "Implicit_Extension_Curve.h"
#include "ImplicitHugoniotCurve.h"
#include "HugoniotContinuation2D2D.h"
#include "HugoniotContinuation_nDnD.h"
#include "ThreePhaseFlowAccumulation.h"
#include "ThreePhaseFlowPermeability.h"
#include "ThreePhaseFlowImplicitHugoniotCurve.h"
//#include "ThreePhaseFlowWaveCurveFactory.h"
#include "WaveCurveFactory.h"
#include "ThreePhaseFlowViscosity.h"
#include "ThreePhaseFlowMobility.h"

class ThreePhaseFlowSubPhysics : public SubPhysics {
    private:
        // No private members.
        //
    protected:
        // This should be moved to ThreePhaseFlowBoundary.
        //
        RealVector W_vertex, O_vertex, G_vertex;

        Implicit_Extension_Curve *iec;

        Parameter *muw_parameter, *muo_parameter, *mug_parameter;
        Parameter *grw_parameter, *gro_parameter, *grg_parameter;
        Parameter *vel_parameter;

        ThreePhaseFlowPermeability *permeability_;
        ThreePhaseFlowViscosity *viscosity_;
        ThreePhaseFlowMobility *mobility_;
    public:
        ThreePhaseFlowSubPhysics();
        virtual ~ThreePhaseFlowSubPhysics();

        // Return some parameters.
        //
        RealVector W(){return W_vertex;}
        RealVector O(){return O_vertex;}
        RealVector G(){return G_vertex;}

//        // These three methods below will, in all probability, move to CoreyQuadSubPhysics.
//        //
//        RealVector E(){
//            RealVector Ep(2);
//            Ep(0) = 0.0;
//            Ep(1) = muo_parameter->value()/(muo_parameter->value() + mug_parameter->value());

//            return Ep;
//        }

//        RealVector B(){
//            RealVector Bp(2);
//            Bp(0) = muw_parameter->value()/(muw_parameter->value() + mug_parameter->value());
//            Bp(1) = 0.0;

//            return Bp;
//        }

//        RealVector D(){
//            RealVector Dp(2);
//            Dp(0) = muw_parameter->value()/(muw_parameter->value() + muo_parameter->value());
//            Dp(1) = muo_parameter->value()/(muw_parameter->value() + muo_parameter->value());

//            return Dp;
//        }

        Parameter* muw(){return muw_parameter;}
        Parameter* muo(){return muo_parameter;}
        Parameter* mug(){return mug_parameter;}
        Parameter* grw(){return grw_parameter;}
        Parameter* gro(){return gro_parameter;}
        Parameter* grg(){return grg_parameter;}
        Parameter* vel(){return vel_parameter;}

        // Permeability.
        //
        ThreePhaseFlowPermeability* permeability(){return permeability_;}

        // Viscosity.
        //
        ThreePhaseFlowViscosity* viscosity(){return viscosity_;}

        // Mobility.
        //
        ThreePhaseFlowMobility* mobility(){return mobility_;}
};

#endif // _THREEPHASEFLOWSUBPHYSICS_


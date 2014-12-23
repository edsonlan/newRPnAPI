#ifndef _THREEPHASEFLOWMOBILITY_
#define _THREEPHASEFLOWMOBILITY_

class ThreePhaseFlowSubPhysics;
class ThreePhaseFlowPermeability;
class ThreePhaseFlowViscosity;

#include "RealVector.h"

class ThreePhaseFlowMobility {
    private:
    protected:
        ThreePhaseFlowPermeability *permeability;
        ThreePhaseFlowViscosity *viscosity;
        ThreePhaseFlowSubPhysics *subphysics;
    public:
        ThreePhaseFlowMobility(ThreePhaseFlowSubPhysics *s);
        virtual ~ThreePhaseFlowMobility();

        virtual double water_mobility(const RealVector &p);
        virtual double oil_mobility(const RealVector &p);
        virtual double gas_mobility(const RealVector &p);
};

#endif // _THREEPHASEFLOWMOBILITY_


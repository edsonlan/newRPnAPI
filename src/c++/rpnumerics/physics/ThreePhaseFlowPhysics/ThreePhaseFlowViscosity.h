#ifndef _THREEPHASEFLOWVISCOSITY_
#define _THREEPHASEFLOWVISCOSITY_

#define VISCOSITY_OK    0
#define VISCOSITY_ERROR 1

#include "WaveState.h"
#include "JetMatrix.h"

class ThreePhaseFlowSubPhysics;

class ThreePhaseFlowViscosity {
    private:
    protected:
        ThreePhaseFlowSubPhysics *subphysics_;
    public:
        ThreePhaseFlowViscosity(ThreePhaseFlowSubPhysics *t);
        virtual ~ThreePhaseFlowViscosity();

        virtual int gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet) = 0;

        int gas_viscosity_jet(const RealVector &w, int degree, JetMatrix &mug_jet){
            return gas_viscosity_jet(WaveState(w), degree, mug_jet);
        }

        // Viscosities.
        //
        virtual double water_viscosity(const RealVector &p);
        virtual double oil_viscosity(const RealVector &p);
        virtual double gas_viscosity(const RealVector &p);
};

#endif // _THREEPHASEFLOWVISCOSITY_


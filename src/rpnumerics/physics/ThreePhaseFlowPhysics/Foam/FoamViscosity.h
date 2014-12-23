#ifndef _FOAMVISCOSITY_
#define _FOAMVISCOSITY_

#include "Parameter.h"
#include "ThreePhaseFlowViscosity.h"

class FoamViscosity: public ThreePhaseFlowViscosity {
    private:
    protected:
        Parameter *mug0_parameter; // mug0

        Parameter *epdry_parameter;
        Parameter *fdry_parameter;
        Parameter *foil_parameter;
        Parameter *fmdry_parameter;
        Parameter *fmmob_parameter;
        Parameter *fmoil_parameter;

        void Fdry(double x, int degree, JetMatrix &fdry_jet);
        void Fo(double so, int degree, JetMatrix &fo_jet);
    public:
        FoamViscosity(Parameter *mug0, 
                      Parameter *epdry,
                      Parameter *fdry,
                      Parameter *foil,
                      Parameter *fmdry,
                      Parameter *fmmob,
                      Parameter *fmoil,
                      ThreePhaseFlowSubPhysics *t);
        ~FoamViscosity();

        int gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet);

        double gas_viscosity(const RealVector &p);
};

#endif // _FOAMVISCOSITY_


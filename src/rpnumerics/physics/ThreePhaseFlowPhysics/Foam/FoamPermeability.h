#ifndef _FOAMPERMEABILITY_
#define _FOAMPERMEABILITY_

#include "ThreePhaseFlowPermeability.h"

class FoamPermeability : public ThreePhaseFlowPermeability {
    private:
    protected:
        Parameter *cnw_parameter, *cno_parameter, *cng_parameter;
        Parameter *nw_parameter, *no_parameter, *ng_parameter;

//        double power(double x, double y);

//        double pow(double x, double y){return 1.0;}
    public:
        FoamPermeability(Parameter *cnw, Parameter *cno, Parameter *cng,
                         Parameter  *nw, Parameter  *no, Parameter  *ng, 
                         ThreePhaseFlowSubPhysics *s);
        virtual ~FoamPermeability();

        int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water, double &reduced_water);
        int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas, double &reduced_gas);
        int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil, double &reduced_oil);

        int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water);
        int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas);
        int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil);

        void reduced_permeability(const RealVector &state, RealVector &rp);
};

#endif // _FOAMPERMEABILITY_


#ifndef _BROOKS_COREYPERMEABILITY_
#define _BROOKS_COREYPERMEABILITY_

#include "ThreePhaseFlowPermeability.h"

class Brooks_CoreyPermeability : public ThreePhaseFlowPermeability {
    private:
    protected:
        Parameter *lambda_parameter;

        Parameter *cnw_parameter, *cno_parameter, *cng_parameter;

        double power(double x, double y);
    public:
        Brooks_CoreyPermeability(Parameter *lambda, Parameter *cnw, Parameter *cno, Parameter *cng, ThreePhaseFlowSubPhysics *s);
        virtual ~Brooks_CoreyPermeability();

        int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water, double &reduced_water);
        int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas, double &reduced_gas);
        int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil, double &reduced_oil);

        int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water);
        int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas);
        int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil);

        void reduced_permeability(const RealVector &state, RealVector &rp);
};

#endif // _BROOKS_COREYPERMEABILITY_


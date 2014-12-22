#ifndef _KOVALFLUXFUNCTION_
#define _KOVALFLUXFUNCTION_

#include "ThreePhaseFlowFluxFunction.h"

class KovalFluxFunction: public ThreePhaseFlowFluxFunction {
    private:
        Parameter *grw_parameter_, *grg_parameter_, *gro_parameter_;
        Parameter *muw_parameter_, *mug_parameter_, *muo_parameter_;
        Parameter *vel_parameter_;
        Parameter *phi_parameter_; // Porosity.
    protected:
    public:
        KovalFluxFunction(Parameter *grw, Parameter *gro, Parameter *grg, 
                          Parameter *muw, Parameter *muo, Parameter *mug,
                          Parameter *vel);

        virtual ~KovalFluxFunction();

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _COREY_QUADRATIC_


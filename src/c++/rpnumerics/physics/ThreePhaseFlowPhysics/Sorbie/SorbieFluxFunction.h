#ifndef _SORBIEFLUXFUNCTION_
#define _SORBIEFLUXFUNCTION_

#include "ThreePhaseFlowFluxFunction.h"
#include "SorbiePermeability.h"

class SorbieFluxFunction : public FluxFunction {
    private:
    protected:
        Parameter *muw_parameter, *muo_parameter, *mug_parameter;
        Parameter *grw_parameter, *gro_parameter, *grg_parameter;
        Parameter *vel_parameter;

        SorbiePermeability *perm;
    public:
        SorbieFluxFunction(Parameter *muw, Parameter *muo, Parameter *mug, 
                           Parameter *grw, Parameter *gro, Parameter *grg,
                           Parameter *vel,
                           SorbiePermeability *p);
        virtual ~SorbieFluxFunction();

        int jet(const WaveState &w, JetMatrix &f, int degree) const;

};

#endif // _SORBIEFLUXFUNCTION_


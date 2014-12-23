#ifndef _FOAMFLUXFUNCTION_
#define _FOAMFLUXFUNCTION_

#include <stdio.h>
#include <stdlib.h>
#include "ThreePhaseFlowFluxFunction.h"
#include "FoamPermeability.h"
#include "FoamViscosity.h"

class FoamFluxFunction : public ThreePhaseFlowFluxFunction {
    private:
        Parameter *grw_, *grg_, *gro_;
        Parameter *muw_, *mug_, *muo_;
        Parameter *vel_;

        FoamPermeability *permeability_;
        FoamViscosity    *foamvisc_;
    protected:
    public:
        FoamFluxFunction(Parameter *grw, Parameter *gro, Parameter *grg, 
                         Parameter *muw, Parameter *muo,
                         Parameter *vel,
                         FoamPermeability *fp,
                         FoamViscosity *fv);

        virtual ~FoamFluxFunction();

        int jet(const RealVector &w, JetMatrix &m, int degree) const;
        int jet(const WaveState &w, JetMatrix &m, int degree) const {
            RealVector p(w.stateSpaceDim());
            for (int i = 0; i < w.stateSpaceDim(); i++) p(i) = w(i);

            return jet(p, m, degree);
        }
};

#endif // _FOAMFLUXFUNCTION_


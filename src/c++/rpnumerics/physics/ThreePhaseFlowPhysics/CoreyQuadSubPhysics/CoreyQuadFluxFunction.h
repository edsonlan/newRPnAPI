#ifndef _COREY_QUADRATIC_
#define _COREY_QUADRATIC_

#include <stdio.h>
#include <stdlib.h>
#include "FluxFunction.h"

//#include "CoreyQuadFluxFunction_Params.h"

class CoreyQuadFluxFunction : public FluxFunction {
    private:
//        double grw, grg, gro;
//        double muw, mug, muo;
//        double vel;
//
////        Permeability *perm;
//        double krw_p, krg_p, kro_p;
//        double cnw, cng, cno;

        Parameter *grw_parameter_, *grg_parameter_, *gro_parameter_;
        Parameter *muw_parameter_, *mug_parameter_, *muo_parameter_;
        Parameter *vel_parameter_;
    protected:
    public:
        CoreyQuadFluxFunction(Parameter *grw, Parameter *gro, Parameter *grg, 
                  Parameter *muw, Parameter *muo, Parameter *mug,
                  Parameter *vel);

//        CoreyQuadFluxFunction(const CoreyQuadFluxFunction_Params &);
//        CoreyQuadFluxFunction(const CoreyQuadFluxFunction &);
//        CoreyQuadFluxFunction * clone() const;

        virtual ~CoreyQuadFluxFunction();

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
};

#endif // _COREY_QUADRATIC_


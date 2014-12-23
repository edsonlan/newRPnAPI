/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) StoneFluxFunction.h
 */

#ifndef _StoneFluxFunction_H
#define _StoneFluxFunction_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */


#include "ThreePhaseFlowFluxFunction.h"
#include "StonePermeability.h"
#include <math.h>
#include "Parameter.h"

class StoneFluxFunction : public ThreePhaseFlowFluxFunction {
    private:
    protected:
        StonePermeability *perm_;

        // Inherited from ThreePhaseFlowFluxFunction: *muw_, *mug_, *muo_;

        Parameter *grw_, *grg_, *gro_, *vel_;
    public:
        StoneFluxFunction(Parameter *grw, Parameter *grg, Parameter *gro, Parameter *muw, Parameter *mug, Parameter *muo, Parameter *vel, StonePermeability *sp);
        virtual ~StoneFluxFunction();
    
        int jet(const WaveState &u, JetMatrix &m, int degree) const;

        StonePermeability* permeability(){return perm_;}

        Parameter* grw(){return grw_;}
        Parameter* gro(){return gro_;} 
        Parameter* grg(){return grg_;}

        Parameter* muw(){return muw_;}
        Parameter* muo(){return muo_;} 
        Parameter* mug(){return mug_;}

        Parameter* vel(){return vel_;}    
};

#endif //! _StoneFluxFunction_H

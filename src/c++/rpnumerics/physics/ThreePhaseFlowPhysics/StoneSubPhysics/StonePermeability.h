/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) StonePermeability.h
 */

#ifndef _StonePermeability_H
#define _StonePermeability_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include <math.h>
#include "StonePermParams.h"
#include "JetMatrix.h"
#include "ThreePhaseFlowPermeability.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class StonePermeability : public ThreePhaseFlowPermeability {
    private:
    protected:
        Parameter *expw_parameter,  *expg_parameter, *expo_parameter;
        Parameter *expow_parameter, *expog_parameter;
        Parameter *cnw_parameter, *cng_parameter, *cno_parameter;
        Parameter *lw_parameter, *lg_parameter;
        Parameter *low_parameter, *log_parameter;
        Parameter *epsl_parameter;

        int kowden_jet(double sow, int degree, JetMatrix &kowj);
        int kogden_jet(double sog, int degree, JetMatrix &kogj);
    public:
        StonePermeability(ThreePhaseFlowSubPhysics *s);
        virtual ~StonePermeability();
    
        void Diff_PermabilityWater(double, double, double, double&, double&, double&, double&, double&, double&);
        void Diff_PermabilityOil(double, double, double, double&, double&, double&, double&, double&, double&);
        void Diff_PermabilityGas(double, double, double, double&, double&, double&, double&, double&, double&);

        int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water);
        int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil);
        int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas);

        void reduced_permeability(const RealVector &state, RealVector &rp);

        Parameter* expw(){return expw_parameter;}
        Parameter* expo(){return expo_parameter;}
        Parameter* expg(){return expg_parameter;}

        Parameter* expow(){return expow_parameter;}
        Parameter* expog(){return expog_parameter;}

        Parameter* cnw(){return cnw_parameter;}
        Parameter* cno(){return cno_parameter;}
        Parameter* cng(){return cng_parameter;}

        Parameter* lw(){return lw_parameter;}
        Parameter* lg(){return lg_parameter;}

        Parameter* low(){return low_parameter;}
        Parameter* log(){return log_parameter;}

        Parameter* epsl(){return epsl_parameter;}
};

#endif //! _StonePermeability_H

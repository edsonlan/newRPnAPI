/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Double_Contact_Function.h
 */

#ifndef _Double_Contact_Function_H
#define _Double_Contact_Function_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "ThreeImplicitFunctions.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class Double_Contact_Function:public ThreeImplicitFunctions {

private:

public:
    
      virtual void curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues *lg, int lfam,
                   const FluxFunction *rf, const AccumulationFunction *ra, GridValues *rg, int rfam,
                   std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve)=0;

};

#endif //! _Double_Contact_Function_H

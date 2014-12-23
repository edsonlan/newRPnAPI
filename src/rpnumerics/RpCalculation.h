/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RpCalculation.h
 */

#ifndef _RpCalculation_H
#define _RpCalculation_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RpSolution.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class RpCalculation {

private:

public:
    
virtual RpSolution & calc ()=0;
virtual RpSolution & recalc()=0;


};

#endif //! _RpCalculation_H

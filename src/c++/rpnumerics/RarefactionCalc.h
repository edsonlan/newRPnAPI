/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RarefactionCalc.h
 */

#ifndef _RarefactionCalc_H
#define _RarefactionCalc_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RpCalculation.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class RarefactionCalc: public RpCalculation{

private:

public:


    RpSolution & calc();
    RpSolution & recalc();
    
};

#endif //! _RarefactionCalc_H

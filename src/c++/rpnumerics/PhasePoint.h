/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) PhasePoint.h
 */

#ifndef _PhasePoint_H
#define _PhasePoint_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RealVector.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class PhasePoint : public RealVector {

public:
	
	PhasePoint(const RealVector & phaseCoods); //TODO Fix implementation

	RealVector operator()(void);

};






#endif //! _PhasePoint_H

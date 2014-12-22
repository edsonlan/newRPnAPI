/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) DefaultAccumulationParams.h
 */

#ifndef _DefaultAccumulationParams_H
#define _DefaultAccumulationParams_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RealVector.h"
#include "AccumulationParams.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class DefaultAccumulationParams : AccumulationParams {
/* 
 * @todo How to create DefaultAccParams before complete physics initalization ?
 */

public:
	DefaultAccumulationParams(void);
	AccumulationParams defaultParams(void);

};

inline AccumulationParams DefaultAccumulationParams::defaultParams(void)
{
	return * new DefaultAccumulationParams();
}

//
//inline DefaultAccumulationParams::DefaultAccumulationParams(void)
//	: AccumulationParams(* new RealVector())
//{
	// JAVA CODE! daniel@impa.br
	//   // can´t call RPNUMERICS methods before complete physics initalization !
	//   super (rpnumerics.RPNUMERICS.physicsID(),new RealVector(rpnumerics.RPNUMERICS.domainDim()));

//}


#endif //! _DefaultAccumulationParams_H

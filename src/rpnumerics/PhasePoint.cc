/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) PhasePoint.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "PhasePoint.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */



PhasePoint::PhasePoint(const RealVector & phaseCoords)
	: RealVector(phaseCoords)
{
}

 RealVector PhasePoint::operator()(void)
{
	return * this;
}

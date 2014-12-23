/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) WavePoint.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "WavePoint.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */



WavePoint::WavePoint(const int size, const double sigma):Point(RealVector(size)),
	  sigma_(sigma)
{
}

WavePoint::WavePoint(const RealVector coords, const double sigma):Point(coords), sigma_(sigma)


{
}

WavePoint::WavePoint(const int size, double * coords, const double sigma): Point(RealVector(size, coords)),
        sigma_(sigma)
{
}

WavePoint::WavePoint(const WavePoint & copy):Point(copy.getCoord()),sigma_(copy.sigma_)
{
}

double WavePoint::sigma(void)
{
	return sigma_;
}

void WavePoint::sigma(double sigma)
{
	sigma_ = sigma;
}

void WavePoint::operator=(double sigma)
{
	sigma_ = sigma;
}

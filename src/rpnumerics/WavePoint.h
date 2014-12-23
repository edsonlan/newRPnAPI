/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) WavePoint.h
 */

#ifndef _WavePoint_H
#define _WavePoint_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RealVector.h"
#include "Point.h"
/*
 * ---------------------------------------------------------------
 * Definitions:
 */

/*!

TODO:
	To implement an exception in 'operator=()' method.
NOTE : 

@ingroup rpnumerics
*/

class WavePoint : public Point{

private:
	double sigma_;

public:
	WavePoint(const int size, const double sigma);
	WavePoint(const RealVector coords, const double sigma);
	WavePoint(const int size, double * coords, const double sigma);
	WavePoint(const WavePoint & copy);

	double sigma(void);
	void sigma(double sigma);
	void operator=(double sigma);
	void operator=(RealVector & coords);
	void operator=(WavePoint & copy);

};



#endif //! _WavePoint_H

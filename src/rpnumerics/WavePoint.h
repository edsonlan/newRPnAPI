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
#include "WavePointInfo.h"
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

    
    WavePointInfo * info_;
public:

        
	WavePoint(const RealVector coords, WavePointInfo * );


};



#endif //! _WavePoint_H

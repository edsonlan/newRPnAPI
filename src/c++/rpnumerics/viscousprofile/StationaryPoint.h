/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) StationaryPoint.h
 */

#ifndef _StationaryPoint_H
#define _StationaryPoint_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "eigen.h"
#include "RealVector.h"
#include <vector>

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class StationaryPoint {
private:

    RealVector * coords_;
    int type_;

public:


    StationaryPoint(const RealVector &, int);
    const RealVector & coords();
    int type();

    virtual ~StationaryPoint();


};

#endif //! _StationaryPoint_H

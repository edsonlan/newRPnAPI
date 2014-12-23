/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) StonePermParams.h
 */

#ifndef _StonePermParams_H
#define _StonePermParams_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */

#include "RealVector.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class StonePermParams {
private:
    RealVector * comp;

public:

    StonePermParams(double expw, double expg, double expo,
            double expow, double expog,
            double cnw, double cng, double cno,
            double lw, double lg,
            double low, double log,
            double epsl);
    StonePermParams();
    StonePermParams(const StonePermParams &);
    virtual ~StonePermParams();

    void reset();

    double component(int);

    const RealVector & params() const;
};

inline const RealVector & StonePermParams::params()const {
    return *comp;
}


#endif //! _StonePermParams_H

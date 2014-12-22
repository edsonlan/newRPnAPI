/* 
 * File:   SegmentedCurve.h
 * Author: edsonlan
 *
 * Created on December 22, 2014, 11:49 AM
 */

#ifndef SEGMENTEDCURVE_H
#define	SEGMENTEDCURVE_H

#include "MyCurve.h"
#include "Segment.h"

class SegmentedCurve: public MyCurve {
public:
    SegmentedCurve(list <Segment *> *);
    SegmentedCurve(const SegmentedCurve& orig);
    virtual ~SegmentedCurve();
private:

};

#endif	/* SEGMENTEDCURVE_H */


/* 
 * File:   PointCurve.h
 * Author: edsonlan
 *
 * Created on December 19, 2014, 12:07 PM
 */

#ifndef POINTCURVE_H
#define	POINTCURVE_H
#include "MyCurve.h"
#include "Point.h"


class PointCurve :public MyCurve {
public:
    PointCurve(list<Point *> *);
    PointCurve(const PointCurve& orig);
    virtual ~PointCurve();
private:

};

#endif	/* POINTCURVE_H */


/* 
 * File:   WavePointInfo.h
 * Author: edsonlan
 *
 * Created on January 5, 2015, 11:27 AM
 */

#ifndef WAVEPOINTINFO_H
#define	WAVEPOINTINFO_H


#include "MyCurve.h"
#include "Point.h"

class WavePointInfo {
    
    double sigma_;
    double * eigenValues_;
    int eigenValuesLength_;
    
    Point * originPoint_;
    MyCurve * originCurve_;
    
    
    
public:
    WavePointInfo(double, int,double *,Point *,MyCurve *);
    virtual ~WavePointInfo();
    
private:

};

#endif	/* WAVEPOINTINFO_H */


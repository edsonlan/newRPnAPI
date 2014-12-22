/* 
 * File:   SegmentedCurveView.h
 * Author: edsonlan
 *
 * Created on September 26, 2014, 2:20 PM
 */

#ifndef SEGMENTEDCURVEVIEW_H
#define	SEGMENTEDCURVEVIEW_H


#include "RealVector.h"
#include "CurveView.h"
#include <vector>

using namespace std;


class SegmentedCurveView:public CurveView {
    
    
public:
    SegmentedCurveView(const vector<vector<RealVector> > &);
    SegmentedCurveView(const SegmentedCurveView& orig);
    void draw (const Graphics *);
    
    virtual ~SegmentedCurveView();
private:
    
    
    vector<vector<RealVector> >  * segments_ ;

};

#endif	/* SEGMENTEDCURVEVIEW_H */


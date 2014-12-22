/* 
 * File:   CurveView.h
 * Author: edsonlan
 *
 * Created on September 26, 2014, 2:25 PM
 */

#ifndef CURVEVIEW_H
#define	CURVEVIEW_H


#include "Graphics.h"
#include"RealVector.h"

class CurveView {
public:
    void virtual draw(Graphics *) = 0;


    virtual ~CurveView();


};

inline CurveView::~CurveView() {

}


#endif	/* CURVEVIEW_H */


/* 
 * File:   PointCurve.cpp
 * Author: edsonlan
 * 
 * Created on December 19, 2014, 12:07 PM
 */

#include "PointCurve.h"

PointCurve::PointCurve(list<Point *> * points):MyCurve((list<Data *> *)points) {
}

PointCurve::PointCurve(const PointCurve& orig):MyCurve(new list<Data *>(*orig.getElements())) {
}



PointCurve::~PointCurve() {

}


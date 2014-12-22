/* 
 * File:   SegmentedCurve.cpp
 * Author: edsonlan
 * 
 * Created on December 22, 2014, 11:49 AM
 */

#include "SegmentedCurve.h"


SegmentedCurve::SegmentedCurve(list <Segment *> * segments):MyCurve((list<Data *> *)segments) {
}

SegmentedCurve::SegmentedCurve(const SegmentedCurve& orig):MyCurve(new list<Data *>(*orig.getElements())) {
}

SegmentedCurve::~SegmentedCurve() {
}


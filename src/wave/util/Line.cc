#include "Line.h"

Line::Line(const RealVector &pp0, const RealVector &pp1): p0(pp0), p1(pp1){
}

Line::~Line(){
}

bool Line::intersect(const Line &l, RealVector &p){
    double alpha, beta;

    return segment_segment_intersection(p0, p1, l.p0, l.p1, p, alpha, beta);
}


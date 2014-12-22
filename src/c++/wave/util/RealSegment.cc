#include "RealSegment.h"

RealSegment::RealSegment(const RealVector & p1, const RealVector & p2) : 
	p1_(new RealVector(p1)), 
	p2_(new RealVector(p2)) 
{
}

RealSegment::RealSegment(const RealSegment & copy) : 
	p1_(copy.p1_), 
	p2_(copy.p2_) 
{
}



RealSegment & RealSegment::operator=(const RealSegment & copy) {

    p1_ = copy.p1_;
    p2_ = copy.p2_;
    
    return *this;
}

bool RealSegment::operator ==(const RealSegment &test) {
    
    if (p1_ == test.p1_ && p2_ == test.p2_) {
        return true;
    }
    else {
        return false;
    }
}

RealSegment::~RealSegment()
{
  delete p1_;
  delete p2_;
}


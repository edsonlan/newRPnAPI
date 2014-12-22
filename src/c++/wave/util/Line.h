#ifndef _LINE_
#define _LINE_

#include "RealVector.h"

class Line {
    private:
    protected:
        RealVector p0, p1;
    public:
        Line(const RealVector &pp0, const RealVector &pp1);
        virtual ~Line();

        bool intersect(const Line &l, RealVector &p);
};

#endif // _LINE_


#ifndef _BOXND_
#define _BOXND_

#include "PointND.h"

class BoxND {
    private:
    protected:
        double max(double x, double y) const {return (x > y) ? x : y;}
        double min(double x, double y) const {return (x < y) ? x : y;}
        double abs(double x) const {return (x > 0.0) ? x : (-x);}
    public:
        PointND pmin, pmax;

        BoxND();
        BoxND(const PointND &ppmin, const PointND &ppmax);
        //BoxND(PointND &ppmin, PointND &ppmax);
        BoxND(const BoxND &b);
        BoxND(const BoxND *b);

        ~BoxND();

        bool intersect(const BoxND &b) const ;
};

#endif // _BOXND_


#ifndef _ZEROIMPLICIT_
#define _ZEROIMPLICIT_

#include <vector>
#include "RealVector.h"

class ZeroImplicit {
    private:
    protected:
    public:
        ZeroImplicit();
        virtual ~ZeroImplicit();

        virtual void find_zeroes(int side, std::vector<RealVector> &zeroes) = 0;
};

#endif // _ZEROIMPLICIT_


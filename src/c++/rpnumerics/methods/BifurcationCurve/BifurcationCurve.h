#ifndef _BIFURCATIONCURVE_
#define _BIFURCATIONCURVE_

#include <vector>
#include <string>
#include "RealVector.h"

class BifurcationCurve {
    private:
    protected:
    public:
        BifurcationCurve();
        virtual ~BifurcationCurve();

        virtual void list_of_secondary_bifurcation_curves(std::vector<int> &type, 
                                                          std::vector<std::string> &name, 
                                                          std::vector<void*> &object,
                                                          std::vector<double (*)(void*, const RealVector &)> &function) = 0;

        virtual void curve(int type, std::vector<RealVector> &c) = 0;

};

#endif // _BIFURCATIONCURVE_


#ifndef _CLASSIFYSHOCKCURVE_
#define _CLASSIFYSHOCKCURVE_

#include <string>
#include "RealVector.h"
#include "ReferencePoint.h"
#include "HugoniotContinuation.h"

class ClassifyShockCurve {
    private:
    protected:
        const HugoniotContinuation *hc;

        std::string sp, sm, sd;

        virtual int half_classify(double sigma, const std::vector<eigenpair> &e, std::string &s);
    public:
        ClassifyShockCurve(const HugoniotContinuation *h);
        virtual ~ClassifyShockCurve();

        int classify_point(const RealVector &p, const ReferencePoint &ref, std::string &s);
};

#endif // _CLASSIFYSHOCKCURVE_


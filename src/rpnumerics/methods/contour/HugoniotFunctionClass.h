#ifndef _HUGONIOTFUNCTIONCLASS_
#define _HUGONIOTFUNCTIONCLASS_

#include "RealVector.h"
#include <vector>
#include "FluxFunction.h"
#include "AccumulationFunction.h"

class HugoniotFunctionClass {
private:
    RealVector * uRef_;
    const FluxFunction *fluxFunction_;

protected:
public:

    HugoniotFunctionClass();
    HugoniotFunctionClass(const FluxFunction &);
    //    HugoniotFunctionClass(FluxFunction *);
    virtual double HugoniotFunction(const RealVector &u) = 0;
    virtual void completeCurve(std::vector<RealVector> &);
    RealVector & getReferenceVector();
    virtual void setReferenceVector(const RealVector &);
    const FluxFunction & getFluxFunction()const;
    virtual void setFluxFunction(const FluxFunction *);
    virtual void setAccumulationFunction(const AccumulationFunction *);

    virtual ~HugoniotFunctionClass();

};

#endif // _HUGONIOTFUNCTIONCLASS_


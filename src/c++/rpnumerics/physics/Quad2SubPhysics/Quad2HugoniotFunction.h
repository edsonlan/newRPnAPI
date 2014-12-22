#ifndef _QUAD2HUGONIOTFUNCTION_
#define _QUAD2HUGONIOTFUNCTION_

#include "Quad2FluxFunction.h"
#include "HugoniotFunctionClass.h"
#include <vector>
#include "eigen.h" // TODO: Find the place

class Quad2HugoniotFunction : public HugoniotFunctionClass {
    private:

        RealVector Uref;
        JetMatrix  UrefJetMatrix;

        std::vector<eigenpair> ve_uref;

        bool Uref_is_elliptic;
    protected:
    public:
        Quad2HugoniotFunction(const RealVector &U, const Quad2FluxFunction &);
        void setReferenceVector(const RealVector & refVec);
        ~Quad2HugoniotFunction();

        double HugoniotFunction(const RealVector &u);
};

#endif // _QUAD2HUGONIOTFUNCTION_


#ifndef _COREY_QUADRATIC_PARAMS_
#define _COREY_QUADRATIC_PARAMS_

#include "FluxParams.h"

class CoreyQuad_Params : public FluxParams {
    private:
    protected:
    public:
        CoreyQuad_Params(const double grw, const double grg, const double gro,
                     const double muw, const double mug, const double muo,
                     const double vel);
        CoreyQuad_Params();
        CoreyQuad_Params(const RealVector &);
        CoreyQuad_Params(const CoreyQuad_Params &);

        ~CoreyQuad_Params();
};

#endif // _COREY_QUADRATIC_PARAMS_


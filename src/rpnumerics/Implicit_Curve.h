#ifndef _IMPLICIT_CURVE_
#define _IMPLICIT_CURVE_

#include "HugoniotFunctionClass.h"
#include "ImplicitFunction.h"

class Implicit_Curve : public ImplicitFunction {
    private:
    protected:
    public:
        Implicit_Curve();
        ~Implicit_Curve();

        int valueFunction(double *foncub, int i, int j, int is_square);
};

#endif // _IMPLICIT_CURVE_


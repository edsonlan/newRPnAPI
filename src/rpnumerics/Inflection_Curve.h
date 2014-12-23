#ifndef _INFLECTION_CURVE_
#define _INFLECTION_CURVE_

#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Matrix.h"
#include "RealVector.h"
#include <vector>
#include "eigen.h"
#include "ImplicitFunction.h"

class Inflection_Curve : public ImplicitFunction {
    private:
    protected:
        int family;
    public:
        Inflection_Curve(){gv = 0;}
        ~Inflection_Curve(){}

        int function_on_square(double *foncub, int i, int j);

        int consistency(double *v1, double *v2, int &orient);

        int curve(const FluxFunction *f, const AccumulationFunction *a, 
                  GridValues &g, int fam, std::vector<RealVector> &inflection_curve);
};

#endif // _INFLECTION_CURVE_


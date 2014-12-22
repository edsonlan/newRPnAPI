#ifndef _ELLIPTIC_BOUNDARY_
#define _ELLIPTIC_BOUNDARY_

#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Matrix.h"
#include "RealVector.h"

#include "ImplicitFunction.h"
#include "ContourMethod.h"
#include "ColorCurve.h"
#include "GridValues.h"

class Elliptic_Boundary : public ImplicitFunction {
    private:
        const FluxFunction *ff;
        const AccumulationFunction *aa;
    protected:
    public:
        Elliptic_Boundary(){gv = 0;}
        ~Elliptic_Boundary();

        int function_on_square(double *foncub, int i, int j);

        int curve(const FluxFunction *f, const AccumulationFunction *a, 
                  GridValues &g, std::vector<RealVector> &elliptic_boundary);

        void map(const RealVector &p, double &f, RealVector &map_Jacobian);

        bool improvable(void);
};

#endif // _ELLIPTIC_BOUNDARY_

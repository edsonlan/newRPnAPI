#ifndef _EQUATIONFUNCTIONLEVELCURVE_
#define _EQUATIONFUNCTIONLEVELCURVE_

#include "ImplicitFunction.h"
#include "ContourMethod.h"

#include <deque>

// TODO: This line below and the ones to it related will change
//       to accomodate the use of Equation (the class, I mean).
//
#include "RpFunction.h"

class EquationFunctionLevelCurve: public ImplicitFunction {
    private:
    protected:
        const RpFunction *function_;

        int component_;

        double level_;

        static double base_level_function(EquationFunctionLevelCurve *obj, const RealVector &p);
        double (*level_function)(EquationFunctionLevelCurve *obj, const RealVector &p);

        virtual double level(const RealVector &p);
    public:
        EquationFunctionLevelCurve(const RpFunction *rpf, GridValues *g);
        virtual ~EquationFunctionLevelCurve();

        virtual int function_on_square(double *foncub, int i, int j);
        virtual void curve(const RealVector &ref, int component, std::vector<RealVector> &c);

        virtual double level(const RealVector &ref, int component);
};

#endif // _EQUATIONFUNCTIONLEVELCURVE_


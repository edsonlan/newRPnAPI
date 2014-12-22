#ifndef _NEWTON_IMPROVEMENT_
#define _NEWTON_IMPROVEMENT_

#ifndef NEWTON_IMPROVEMENT_OK
#define NEWTON_IMPROVEMENT_OK 0
#endif

#ifndef NEWTON_IMPROVEMENT_ERROR
#define NEWTON_IMPROVEMENT_ERROR 1
#endif

#include "RealVector.h"
#include "Matrix.h"
#include "ImplicitFunction.h"
#include <stdio.h>

class Newton_Improvement {
    private:
    protected:
        const ImplicitFunction *imp_map;

        double alpha_init(const RealVector &p0, const RealVector &p1, const RealVector &p);
        void function_and_derivative(double alpha, const RealVector &p0, const RealVector &p1, double &f_alpha, double &df_alpha);
    public:
        Newton_Improvement(const ImplicitFunction *m); 
        ~Newton_Improvement();
        int newton(const RealVector &p0, const RealVector &p1, const RealVector &p_init, RealVector &p_improved);
};

#endif // _NEWTON_IMPROVEMENT_


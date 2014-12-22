#ifndef _EULERSOLVER_
#define _EULERSOLVER_

#include "ODE_Solver.h"
#include "RealVector.h"
#include "Bisection.h"
#include "Boundary.h"

#include <algorithm> // For std::max()

class EulerSolver : public ODE_Solver {
    private:
    protected:
        const Boundary *boundary;

        int number_of_subdivisions;
    public:
        EulerSolver(const Boundary *b, int n);
        virtual ~EulerSolver();

        int integrate_step(int (*field)(int*, double*, double*, double*, int*, double*),
                           int (*jacobianfield)(int *, double *, double *, int *, int *, double *, int *),  
                           int *function_object, double *function_data, 
                           const double init_time,  const RealVector &init_point,  
                           const double final_time,       RealVector &final_point/*,
                           int *es_object, int *es_data*/) const;
};

#endif // _EULERSOLVER_


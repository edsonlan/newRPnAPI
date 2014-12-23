#ifndef _LSODE_
#define _LSODE_

#include "ODE_Solver.h"

extern "C"{
    int lsode_(int (*)(int *, double *, double *, double *, int *, double *), int *, double *, double *, double *,
               int *, double *, double *, int *, int *, int *, double *, int *,
               int *, int *, int(*)(int *, double *, double *, int *, int *, double *, int *), int *, int*, double*);
}

class LSODE : public ODE_Solver {
    private:
    protected:
        //int istate;
    public:
        LSODE();
        virtual ~LSODE();

        int integrate_step(int (*field)(int*, double*, double*, double*, int*, double*),
                           int (*jacobianfield)(int *, double *, double *, int *, int *, double *, int *),  
                           int *function_object, double *function_data, 
                           const double init_time,  const RealVector &init_point,  
                           const double final_time,       RealVector &final_point/*,
                           int *lsode_object, int *lsode_data*/) const;

        
};

#endif // _LSODE_


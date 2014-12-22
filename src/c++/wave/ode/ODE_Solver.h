#ifndef _ODE_SOLVER_
#define _ODE_SOLVER_

#include "ODE_Solver_Return_Codes.h"
#include "RealVector.h"

// Base class for ODE solvers.
//
class ODE_Solver {
    private:
    protected:
    public:
        virtual ~ODE_Solver(){}

        virtual int integrate_step(int (*field)(int*, double*, double*, double*, int*, double*), int (*jacobian_field)(int *, double *, double *, int *, int *, double *, int *),
                                   int *function_object, double *function_data, 
                                   const double init_time,  const RealVector &init_point,  
                                   const double final_time,       RealVector &final_point/*,
                                   int *ode_object, int *ode_data*/) const = 0;

        // This method can be used if the Jacobian of the field is unknown.
        //
        virtual int integrate_step(int (*field)(int*, double*, double*, double*, int*, double*), 
                                   int *function_object, double *function_data,
                                   const double init_time,  const RealVector &init_point,  
                                   const double final_time,       RealVector &final_point/*,
                                   int *ode_object, int *ode_data*/) const {

            return integrate_step(field, 0, 
                                  function_object, function_data, 
                                  init_time,  init_point,  
                                  final_time, final_point/*,
                                  ode_object, ode_data*/);
        }

        // This method can be used to initialize the solver.
        // For LSODE in particular there is a difference between the first iteration and the rest.
        //
        virtual void initialize(){
            return;
        }
};

#endif // _ODE_SOLVER_


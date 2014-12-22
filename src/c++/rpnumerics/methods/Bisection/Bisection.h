#ifndef _BISECTION_
#define _BISECTION_

#include <cmath>
#include "RealVector.h"
#include "ODE_Solver.h"

#define BISECTION_CONVERGENCE_OK    0
#define BISECTION_CONVERGENCE_ERROR 1

#define BISECTION_FUNCTION_OK       0
#define BISECTION_FUNCTION_ERROR    2

#define BISECTION_SOLVER_OK         0
#define BISECTION_SOLVER_ERROR      3

#define BISECTION_EQUAL_SIGN        4

// Given a function f, and an interval [a, b] (or [b, a]), such that f(a)*f(b) < 0.0, this
// method finds a value c in the interval, such that f(c) = 0.0, approximately,
// by means of the bisection method.
//
//template <class double, class type_out>
class Bisection {
    private:
    protected:
    public:
        static int bisection_method(const double &t_in,  const RealVector &p_in,
                                const double &t_fin, const RealVector &p_fin,
                                double epsilon, 
                                double &c_t, RealVector &p_c,
                                //int (*f)(const RealVector &x, RealVector &y, int *f_o, int *f_d), int *function_object, int *function_data,
                                int (*f)(int*, double*, double*, double*, int*, double*), int *function_object, double *function_data,
                                //int (*ode_solver)(int (*field)(const RealVector &, RealVector &, int *, int *), int * /*function_object*/, int * /*function_data*/,
                                const ODE_Solver *odesolver,
                                //int (*ode_solver)(int (*field)(int*, double*, double*, double*, int*, double*), int * /*function_object*/, double * /*function_data*/,  
                                //                  const double init_time,  const RealVector &init_point,  
                                //                  const double final_time,       RealVector &final_point,
                                //                  int *, int *), int *ode_solver_object, int *ode_solver_data,
                                int (*signal_event)(const RealVector & where, double & event_measure, int *seo, int *sed), int *signal_event_object, int *signal_event_data);

        static int simple_bisection(double x_init, double x_end, void *obj, int (*)(void*, double, double&), double &x0);
};

#endif // _Bisection_


#ifndef _RAREFACTION_
#define _RAREFACTION_

#include "RarefactionCurve.h"

//#ifndef
//#define RAREFACTION_INIT_OK        0
//#endif

#define RAREFACTION_INIT_FAILURE   1

#define RAREFACTION_NOT_MONOTONOUS 2

#define RAREFACTION_SIMPLE_ACCUMULATION  10

#define RAREFACTION_GENERAL_ACCUMULATION 11

//#ifndef
//#define RAREFACTION_SPEED_INCREASE 20
//#endif

//#ifndef
//#define RAREFACTION_SPEED_NEUTRAL  21
//#endif

//#ifndef
//#define RAREFACTION_SPEED_DECREASE 22
//#endif

//#ifndef
//#define RAREFACTION_INITIALIZE_YES 30
//#endif

//#ifndef
//#define RAREFACTION_INITIALIZE_NO  31
//#endif

#define RAREFACTION_FOR_ITSELF                   40

#define RAREFACTION_AS_ENGINE_FOR_INTEGRAL_CURVE 41

//#ifndef
//#define RAREFACTION_REACHED_BOUNDARY             100
//#endif

//#ifndef
//#define RAREFACTION_REACHED_INFLECTION           101
//#endif

#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "eigen.h"
#include "Boundary.h"
#include <vector>
#include <fstream>

#include "ODE_Solver.h"
#include "Bisection.h"
#include "EulerSolver.h"

//#define TEST_STENCIL

extern "C"{
    int lsode_(int (*)(int *, double *, double *, double *, int *, double *), int *, double *, double *, double *,
            int *, double *, double *, int *, int *, int *, double *, int *,
            int *, int *, int(*)(int *, double *, double *, int *, int *, double *, int *), int *, int*, double*);
}

class Rarefaction {
    private:
        static FluxFunction         *fluxfunction;
        static AccumulationFunction *accumulationfunction;
        static const Boundary       *boundary;
        static int type;
        static int family;

        // For test purposes only: delete afterwards, or comment the #define TEST_STENCIL above.
        #ifdef TEST_STENCIL
        static std::ofstream stencil;
        #endif

        static double ddot(int n, double *x, double *y);
        static void matrixmult(int m, int p, int n, double *A, double *B, double *C);

        static void fill_with_jet(const RpFunction *flux_object, int n, double *in, int degree, double *F, double *J, double *H);

        static int flux(int *neq, double *xi, double *in, double *out, int *nparam, double *param);

        static int compute_last_point(const RealVector &previous_point, const RealVector &new_point, RealVector &last_point);

        static void compute_all_eigenpairs(int n, const RealVector &in, std::vector<eigenpair> &e);

        static void compute_eigenpair(int n, const RealVector &in, double &lambda, RealVector &eigenvector);

        static double compute_lambda(int n, const RealVector &in);

        static int init(const RealVector &initial_point, int increase, double deltaxi, RealVector &second_point);

        static double dirdrv(int n, const RealVector &p, const RealVector &direction);
        static double dirdrv(int n, double *point, double *dir);
        static int initial_dirdrv(int n, const RealVector &p, int increase, double &dd, RealVector &direction);

        static int rar_last_point(int n, const RealVector &p0, const RealVector &p1, RealVector &out);
    protected:
    public:
        static int curve(const RealVector &initial_point, 
                         int initialize,
                         const RealVector *initial_direction,
                         int curve_family, 
                         int increase,
                         int type_of_rarefaction,
                         ODE_Solver *odesolver,
                         double deltaxi,
                         const FluxFunction *ff, const AccumulationFunction *aa,
                         int type_of_accumulation,
                         const Boundary *boundary,
                         std::vector<RealVector> &rarcurve,
                         std::vector<RealVector> &inflection_points);

        static int flux_wrapper(const RealVector &in, RealVector &out, int *flux_object, int *flux_data);

//        static int speed_increase(void){return RAREFACTION_SPEED_INCREASE;}
//        static int speed_neutral(void){return RAREFACTION_SPEED_NEUTRAL;}
//        static int speed_decrease(void){return RAREFACTION_SPEED_DECREASE;}
        static int inflection_signal_event(const RealVector & where, double & directional_derivative_measure, int *signal_object, int *reference_direction);

//    static int lsode_solver(int (*field)(const RealVector &, RealVector &, int *, int *), int *function_object, int *function_data, 
//                              const double init_time,  const RealVector &init_point,  
//                              const double final_time,       RealVector &final_point,
//                              int *lsode_object, int *lsode_data);
};

#endif // _RAREFACTION_


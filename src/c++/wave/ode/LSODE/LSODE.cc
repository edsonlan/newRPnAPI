#include "LSODE.h"

LSODE::LSODE(){
}

LSODE::~LSODE(){
}

int LSODE::integrate_step(int (*field)(int*, double*, double*, double*, int*, double*),
                          int (*jacobianfield)(int *, double *, double *, int *, int *, double *, int *),
                          int *function_object, double *function_data, 
                          const double init_time,  const RealVector &init_point,  
                          const double final_time,       RealVector &final_point/*,
                          int *lsode_object, int *lsode_data*/) const {

    // Dimension of the space:
    int n = init_point.size();

    int ml; // Not used.
    int mu; // Not used.

    // ???
    int nrpd = 4;

    // Is the tolerance the same for all the elements of U (1) or not (2)?
    int itol = 2; // 1: atol scalar; 2: atol array.
    double rtol = 1e-4;
    double *atol = new double[n];
    for (int i = 0; i < n; i++) atol[i] = 1e-6;

    int mf;
    if (jacobianfield != 0) mf = 21; // Jacobian provided by the user.
    else                    mf = 22; // Jacobian NOT provided by the user.

    // Lsode uses rwork to perform its computations.
    // lrw is the declared length of rwork
    int lrw;
    if (mf == 10) lrw = 20 + 16 * n;
    else if (mf == 21 || mf == 22) lrw = 22 + 9*n + n*n;
    else if (mf == 24 || mf == 25) lrw = 22 + 10*n + (2*ml + mu)*n;

    double *rwork = new double[lrw];

    // Normal computation of values at tout.
    int itask = 1;

    // Set to 1 initially.
    // This is where LSODE's info parameter. Must be set to 1 the first time.
    //
    // Since all ODE_Solvers will perform only one step each time (by means of this method), 
    // it would be more difficult to conserve or keep data from one step to the next. Also
    // Bisection will use these ODE_Solvers in a different context, potentially immediatly after
    // computing a point in the curve.
    //
    int istate = 1;

    // No optional inputs
    int iopt = 0;

    // Lsode uses iwork to perform its computations.direction
    // liw is the declared length of iwork
    int liw;
    if (mf == 10) liw = 20;
    else if (mf == 21 || mf == 22 || mf == 24 || mf == 25) liw = 20 + n;

    int *iwork = new int[liw];
    
    // Input and output are in the same place here.
    final_point = init_point;

    double xi     = init_time;
    double new_xi = final_time;

    lsode_(field, &n, final_point.components(), 
           &xi, &new_xi, 
           &itol, &rtol, atol, &itask, 
           &istate, 
           &iopt, 
           rwork, &lrw, 
           iwork, &liw, 
           jacobianfield, 
           &mf, 
           function_object, function_data);

    delete [] iwork;
    delete [] rwork;
    delete [] atol;
    
    // TODO: Perhaps the return codes should be augmented to cover all of LSODE's possible error codes.
    //
    if (istate == ODE_SOLVER_OK)                        return ODE_SOLVER_OK;
    else if (istate == ODE_SOLVER_POINT_OUTSIDE_DOMAIN) return ODE_SOLVER_POINT_OUTSIDE_DOMAIN;
    else                                                return ODE_SOLVER_ERROR;
}


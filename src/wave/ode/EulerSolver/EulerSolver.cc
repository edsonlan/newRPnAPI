#include "EulerSolver.h"

EulerSolver::EulerSolver(const Boundary *b, int n) : boundary(b), number_of_subdivisions(n) {
    number_of_subdivisions = std::max(2, n);
}

EulerSolver::~EulerSolver(){
}

int EulerSolver::integrate_step(int (*field)(int*, double*, double*, double*, int*, double*),
                                int (*jacobianfield)(int *, double *, double *, int *, int *, double *, int *),  
                                int *function_object, double *function_data, 
                                const double init_time,  const RealVector &init_point,  
                                const double final_time,       RealVector &final_point/*,
                                int *es_object, int *es_data*/) const{

    int dimension = init_point.size();

    // Time-related
    double time = init_time;
    double delta_t = (final_time - init_time)/(double)(number_of_subdivisions - 1);

    // Output
    final_point = init_point;

    RealVector out(dimension);
    for (int i = 0; i < number_of_subdivisions; i++){
        if (!boundary->inside(final_point)) return ODE_SOLVER_POINT_OUTSIDE_DOMAIN;

        int info_field = (*field)(&dimension, &time, final_point.components(), out.components(), function_object, function_data);
        if (info_field == FIELD_ERROR) return ODE_SOLVER_ERROR;
        else if (info_field == FIELD_POINT_OUTSIDE_DOMAIN) return ODE_SOLVER_POINT_OUTSIDE_DOMAIN;

        time += delta_t;
        final_point = final_point + delta_t*out;
    }

    return ODE_SOLVER_OK;
}


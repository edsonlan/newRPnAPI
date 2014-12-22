#include <math.h>
#include <cmath>

#include "Newton_Improvement.h"

Newton_Improvement::Newton_Improvement(const ImplicitFunction *m){
    imp_map = m;

    //printf("imp_map = %p\n", imp_map);
}

Newton_Improvement::~Newton_Improvement(){
}

double Newton_Improvement::alpha_init(const RealVector &p0, const RealVector &p1, const RealVector &p){
    int n = p0.size();

    double res = 0.0;
    double segment_squared = 0.0;

    RealVector segment(n);
    for (int i = 0; i < n; i++){
        double segment_component_i = p1.component(i) - p0.component(i);
        segment_squared += segment_component_i*segment_component_i;
        res += segment_component_i*(p.component(i) - p0.component(i));
    }

    return res/segment_squared;
}

void Newton_Improvement::function_and_derivative(double alpha, const RealVector &p0, const RealVector &p1, double &f_alpha, double &df_alpha){
    int n = 2;
    RealVector p(n);
    double beta = 1.0 - alpha;
    for (int i = 0; i < n; i++) p.component(i) = alpha*p1.component(i) + beta*p0.component(i);

    //RealVector map_function(n - 1);
    RealVector map_Jacobian(2);  // n = 2;
    ((ImplicitFunction *)imp_map)->map(p, f_alpha, map_Jacobian);

    // Segment
    RealVector segment(n);
    for (int i = 0; i < n; i++) segment.component(i) = p1.component(i) - p0.component(i);

//    // Function
//    f_alpha = 0.0;
//    for (int i = 0; i < n; i++) f_alpha += map_function(i)*segment.component(i);

    // Derivative
    df_alpha = 0.0;
    for (int i = 0; i < n; i++){
        df_alpha += map_Jacobian(i)*segment.component(i);
    }
        
    return;
} 

int Newton_Improvement::newton(const RealVector &p0, const RealVector &p1, const RealVector &p_init, RealVector &p_improved){
    int max_it = 10; // Maximum number of iterations


    int it = 0;
    double alpha;

    // Segment
//    RealVector segment(n);
//    double segment_squared = 0;
//    for (int i = 0; i < n; i++){
//        segment.component(i) = p1.component(i) - p0.component(i);
//        segment_squared += segment.component(i);
//    }

    alpha = alpha_init(p0, p1, p_init);

    double delta_alpha;

    do {
        double f, df;
        function_and_derivative(alpha, p0, p1, f, df);
        // TODO: The tolerance must be relative, not absolute as it is now.
        if (fabs(df) < 1.0e-10 ) delta_alpha = 0.0;
        else delta_alpha = -f/df;
        alpha += delta_alpha;

        it++;
    } while (it < max_it && fabs(delta_alpha) > 1e-6); // TODO: Find a good value for fabs(delta_alpha), 1e-3 is not necessarily the one.

    p_improved.resize(p0.size());

    if (alpha >= 0.0 && alpha <= 1.0){
        double beta = 1.0 - alpha;
        for (int i = 0; i < p0.size(); i++) p_improved.component(i) = alpha*p1.component(i) + beta*p0.component(i);

        return NEWTON_IMPROVEMENT_OK;
    }
    else {
       for (int i = 0; i < p0.size(); i++) p_improved.component(i) = p_init.component(i);

       printf("Newton did not converge :p alpha = %lf\n", alpha);
       printf("Error = %4.20f, ErrorSQR = %4.20f\n", delta_alpha, delta_alpha*delta_alpha);
       printf("p0 = (%lf, %lf), p1 = (%lf, %lf), p_init = (%lf, %lf)\n", p0.component(0), p0.component(1), p1.component(0), p1.component(1), p_improved.component(0), p_improved.component(1));

       return NEWTON_IMPROVEMENT_ERROR;
    }
}


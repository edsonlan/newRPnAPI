#include "Secondary_Bifurcation.h"

/**
 * The basic idea of the code written here is to exchange the Implicit Function Theorem for an
 * algebraic procedure, see e.g.
 *    W. FULTON (1969) "Algebraic Curves". New York, Benjamin.
 *
 * Thus we can find the auto-intersection of a curve f(x,y) = 0 with the implicit form itself and
 * the two auxiliary constraints:
 *    df(x,y)/dx = 0     and     df(x,y)/dy = 0.
 * In this case f(x,y) takes the RH-condition without speed, i.e., f(x,y) = [f_1][g_2] - [f_2][g_1],
 * where f_i are the flux functions and g_i the accumulation functions. (As always [f] means the
 * difference f(X) - f(X_ref).)
 *
 * Notice that it is mandatory to say that the curve is "singular" because the primary bifurction
 * also satisfies the same loci.
**/

bool Secondary_Bifurcation::function_on_cell(double *val, int ir, int jr, int kl, int kr){
    int domain_i, domain_j;

    if      (kr == 0) {domain_i = ir;     domain_j = jr;}
    else if (kr == 1) {domain_i = ir + 1; domain_j = jr;}
    else if (kr == 2) {domain_i = ir + 1; domain_j = jr + 1;}
    else if (kr == 3) {domain_i = ir;     domain_j = jr + 1;}

    double f1 = gv_right->F_on_grid(domain_i, domain_j).component(0);
    double f2 = gv_right->F_on_grid(domain_i, domain_j).component(1);
    double g1 = gv_right->G_on_grid(domain_i, domain_j).component(0);
    double g2 = gv_right->G_on_grid(domain_i, domain_j).component(1);

    double df11 = gv_right->JF_on_grid(domain_i, domain_j)(0,0);
    double df12 = gv_right->JF_on_grid(domain_i, domain_j)(0,1);
    double df21 = gv_right->JF_on_grid(domain_i, domain_j)(1,0);
    double df22 = gv_right->JF_on_grid(domain_i, domain_j)(1,1);

    double dg11 = gv_right->JG_on_grid(domain_i, domain_j)(0,0);
    double dg12 = gv_right->JG_on_grid(domain_i, domain_j)(0,1);
    double dg21 = gv_right->JG_on_grid(domain_i, domain_j)(1,0);
    double dg22 = gv_right->JG_on_grid(domain_i, domain_j)(1,1);

    double Df1 = f1 - flux_left(0, kl);
    double Df2 = f2 - flux_left(1, kl);
    double Dg1 = g1 - accum_left(0, kl);
    double Dg2 = g2 - accum_left(1, kl);

    // Output        
    val[0] = Df1 * Dg2 - Df2 * Dg1;
    val[1] = df11 * Dg2 + Df1 * dg21 - df21 * Dg1 - Df2 * dg11;
    val[2] = df12 * Dg2 + Df1 * dg22 - df22 * Dg1 - Df2 * dg12;
    
    return true;
}

// ALWAYS: flux and accum are 2x4 matrices.
//
//     3     2
//     o-----o
//     |     |
//     |     |
//     o-----o
//     0     1 
//
//     0 = (i, j),
//     1 = (i + 1, j),
//     2 = (i + 1, j + 1),
//     3 = (i, j + 1).
//
bool Secondary_Bifurcation::prepare_cell(int i, int j) {
    int domain_i, domain_j;

    for (int kr = 0; kr < 4; kr++){
        if      (kr == 0) {domain_i = i;     domain_j = j;}
        else if (kr == 1) {domain_i = i + 1; domain_j = j;}
        else if (kr == 2) {domain_i = i + 1; domain_j = j + 1;}
        else if (kr == 3) {domain_i = i;     domain_j = j + 1;}

        flux_left(0, kr)  = gv_left->F_on_grid(domain_i, domain_j).component(0);
        flux_left(1, kr)  = gv_left->F_on_grid(domain_i, domain_j).component(1);

        accum_left(0, kr) = gv_left->G_on_grid(domain_i, domain_j).component(0);
        accum_left(1, kr) = gv_left->G_on_grid(domain_i, domain_j).component(1);
    }

    return true;
}

void Secondary_Bifurcation::curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues &lg,
                           const FluxFunction *rf, const AccumulationFunction *ra, GridValues &rg,
                           std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve){
    lff = lf;
    laa = la;
    gv_left = &lg;

    rff = rf;
    raa = ra;
    gv_right = &rg;

    gv_left->fill_Jacobians_on_grid(lff, laa);
    gv_right->fill_Jacobians_on_grid(rff, raa);

    left_curve.clear(); 
    right_curve.clear(); 

    Contour2x2_Method::curve2x2(this, left_curve, right_curve);

    return;
}


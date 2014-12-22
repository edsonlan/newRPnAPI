#include "Double_Contact.h"

bool Double_Contact::function_on_cell(double *val, int ir, int jr, int kl, int kr){
    int domain_i, domain_j;

    if ( !(gv_left->cell_is_real(ir, jr)) ) return false;

    if      (kr == 0) {domain_i = ir;     domain_j = jr;}
    else if (kr == 1) {domain_i = ir + 1; domain_j = jr;}
    else if (kr == 2) {domain_i = ir + 1; domain_j = jr + 1;}
    else if (kr == 3) {domain_i = ir;     domain_j = jr + 1;}

//    if (!gv_right->eig_is_real(domain_i, domain_j)[right_family]) return false;

    double lr  = gv_right->e(domain_i, domain_j)[right_family].r;
    double fr  = gv_right->F_on_grid(domain_i, domain_j).component(0);

    double hur = gv_right->G_on_grid(domain_i, domain_j).component(0);
    double gr  = gv_right->F_on_grid(domain_i, domain_j).component(1);
    double hvr = gv_right->G_on_grid(domain_i, domain_j).component(1);

    // Output        
    val[0] = lambda_left[kl] - lr;
    val[1] = lambda_left[kl]*(accum_left(0, kl) - hur) - (flux_left(0, kl) - fr);
    val[2] = lr*(accum_left(1, kl) - hvr) - (flux_left(1, kl) - gr);
    
    return true;
}

// Originally this function was: preplftc, and is defined in: locimp.F.
//
// ALWAYS: flux and accum are 2x4 matrices.
//         lambda is a vector with 4 elements.
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
bool Double_Contact::prepare_cell(int i, int j) {
    int domain_i, domain_j;

    for (int kr = 0; kr < 4; kr++){
        if      (kr == 0) {domain_i = i;     domain_j = j;}
        else if (kr == 1) {domain_i = i + 1; domain_j = j;}
        else if (kr == 2) {domain_i = i + 1; domain_j = j + 1;}
        else if (kr == 3) {domain_i = i;     domain_j = j + 1;}

        if ( !gv_left->eig_is_real(domain_i, domain_j)[left_family] ) return false;
        lambda_left[kr]   = gv_left->e(domain_i, domain_j)[left_family].r;

        flux_left(0, kr)  = gv_left->F_on_grid(domain_i, domain_j).component(0);
        flux_left(1, kr)  = gv_left->F_on_grid(domain_i, domain_j).component(1);

        accum_left(0, kr) = gv_left->G_on_grid(domain_i, domain_j).component(0);
        accum_left(1, kr) = gv_left->G_on_grid(domain_i, domain_j).component(1);
    }

    return true;
}

void Double_Contact::curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues *lg, int lfam,
                           const FluxFunction *rf, const AccumulationFunction *ra, GridValues *rg, int rfam,
                           std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve){
    lff = lf;
    laa = la;
    gv_left = lg;
    left_family = lfam;

    rff = rf;
    raa = ra;
    gv_right = rg;
    right_family = rfam;

    // It is assumed that the grid_value must be the same
    singular = ( (left_family == right_family) && (lg == rg) );

    gv_left->fill_eigenpairs_on_grid(lff, laa);
    gv_right->fill_eigenpairs_on_grid(rff, raa);

    left_curve.clear(); 
    right_curve.clear(); 

    std::cout << "Before" << std::endl;

    Contour2x2_Method::curve2x2(this, left_curve, right_curve);

    std::cout << "After" << std::endl;

    return;
}


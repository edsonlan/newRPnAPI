#include "CompositeCurve.h"
#include "WaveCurve.h"

int CompositeCurve::composite_field(int *two_n, double *xi, double *pointpair, double *field, int *obj, double* /* Not used */){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    // Extract some information about the composite proper:
    //
    int family = composite_object->family;

    const FluxFunction         *f = composite_object->flux;
    const AccumulationFunction *g = composite_object->accum;
    const Boundary             *b = composite_object->boundary;

    // Extract some information about the rarefaction:
    //
    const FluxFunction         *rf = composite_object->rarflux;
    const AccumulationFunction *rg = composite_object->raraccum;
    const Boundary             *rb = composite_object->rarboundary;

    int n = (*two_n)/2;

    // Extract the points and check if they are within the boundary.
    //
    RealVector rarefaction_point(n);
    RealVector composite_point(n);

    for (int i = 0; i < n; i++){
        rarefaction_point(i) = pointpair[i];
        composite_point(i)   = pointpair[i + n];
    }

    if (!rb->inside(rarefaction_point)){
        std::cout << "CompositeCurve::composite_field(): rarefaction point " << rarefaction_point << " is outside the boundary." << std::endl;

        return FIELD_POINT_OUTSIDE_DOMAIN;
    }
    if (!b->inside(composite_point)){
        std::cout << "CompositeCurve::composite_field(): composite point " << composite_point << " is outside the boundary." << std::endl;

        return FIELD_POINT_OUTSIDE_DOMAIN;
    }

    // First compute the rarefaction part of the field: dU^-/dxi = r(U^-) for the given family.
    // Notice that the flux and accumulation used are those associated with the
    // rarefaction.
    //
    RealVector reference_vector = composite_object->reference_vector;

    JetMatrix Fm_jet(n);
    rf->jet(rarefaction_point, Fm_jet, 2);
    DoubleMatrix FmJac = Fm_jet.Jacobian();

    JetMatrix Gm_jet(n);
    rg->jet(rarefaction_point, Gm_jet, 2);
    DoubleMatrix GmJac = Gm_jet.Jacobian();

    std::vector<eigenpair> e;
    
    // TODO: All this info may be stored in the rarefaction, since it had to be computed there anyway.
    Eigen::eig(n, FmJac.data(), GmJac.data(), e);

    if (family < 0 || family > e.size() - 1) return FIELD_ERROR;

    // TODO: Check also that lambda is not complex!

    double lambdam = e[family].r;
    RealVector rm(n, e[family].vrr.data()); // Right eigenvector
    RealVector lm(n, e[family].vlr.data()); // Left  eigenvector

    // Verify that rm points in the right direction.
    // Since rm will be used below, it is necessary that this change is permanent.
    //
    if (rm*reference_vector < 0.0) rm = -rm;

    RealVector Gm = Gm_jet.function();

    // For lambda's directional derivative.
    //
    std::vector<DoubleMatrix> Fm_H = Fm_jet.Hessian();
    std::vector<DoubleMatrix> Gm_H = Gm_jet.Hessian();

    RealVector h(Fm_H.size());

    for (int i = 0; i < Fm_H.size(); i++) h(i) = rm*((Fm_H[i] - lambdam*Gm_H[i])*rm);
    
    double dirdrv = (lm*h)/(lm*(GmJac*rm));

    // Now compute the composite part of the field.
    // Notice that the flux and accumulation used are those associated with the
    // composite.
    //
    JetMatrix Fp_jet(n);
    f->jet(composite_point, Fp_jet, 1);
    
    JetMatrix Gp_jet(n);
    g->jet(composite_point, Gp_jet, 1);
    RealVector Gp = Gp_jet.function();

    DoubleMatrix characteristic_matrix = Fp_jet.Jacobian() - lambdam*Gp_jet.Jacobian();

    if (std::abs(det(characteristic_matrix)) < composite_object->tolerance){
        std::cout << "Composite field. det = " << det(characteristic_matrix) << std::endl;

        return FIELD_ERROR; // TODO: Check this tolerance!
    }

    //
    RealVector cf;
    int info_solve = solve(characteristic_matrix, dirdrv*(Gp - Gm), cf);
    if (info_solve == REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR){
        std::cout << "Composite flux. info_solve = " << info_solve;
        return FIELD_ERROR;
    }

    if (composite_object->normalize_with_respect_to_whom == NORMALIZE_WITH_RESPECT_TO_COMPOSITE){
        double inv_nrm = 1.0/norm(cf);

        rm *= inv_nrm;
        cf *= inv_nrm;
    }

    // The first part of the field:
    //
    for (int i = 0; i < n; i++) field[i] = rm(i);

    // The last part of the field:
     //
    for (int i = 0; i < n; i++) field[i + n] = cf(i);

//    // TODO: This normalization does not appear to improve the quality of the solutions.
//    //       The curves are longer, though, which produces better-looking curves.
//    //
//    //for (int i = 0; i < (*two_n); i++) field[i]first_backwards_rarefaction_point /= nrm;

//    // TODO: This normalization seems to be working ok in terms of results. Nevertheless, since
//    //       the distance between two consecutive points in the rarefaction thus computed is smaller
//    //       than in the original curve, we will not reach the starting point (because).
//    //        
//    // std::cout << "*** Composite. Field = " << RealVector(*two_n, field) << std::endl;

//    // TODO: The last part of the field, corresponding to the composite, MAY need to be normalized.
//    //       In that case the first part should be renormalized.

//    // If the field vector is large, switch to the field near the double contact
//    if (nrm*composite_object->cmp_deltaxi > .005){
//        composite_object->use_field_near_double_contact = true;
//    }

    return FIELD_OK;
}

int CompositeCurve::composite_field_near_double_contact(int *two_n, double *xi, double *pointpair, double *field, int *obj, double* /* Not used */){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    // Extract some information about the composite proper:
    //
    int family = composite_object->family;

    const FluxFunction         *f = composite_object->flux;
    const AccumulationFunction *g = composite_object->accum;
    const Boundary             *b = composite_object->boundary;

    // Extract some information about the rarefaction:
    //
    const FluxFunction         *rf = composite_object->rarflux;
    const AccumulationFunction *rg = composite_object->raraccum;
    const Boundary             *rb = composite_object->rarboundary;

    int n = (*two_n)/2;

    // Extract the points and check if they are within the boundary.
    //
    RealVector rarefaction_point(n);
    RealVector composite_point(n);

    for (int i = 0; i < n; i++){
        rarefaction_point(i) = pointpair[i];
        composite_point(i)   = pointpair[i + n];
    }

    if (!rb->inside(rarefaction_point)){
        std::cout << "CompositeCurve::composite_field(): rarefaction point " << rarefaction_point << " is outside the boundary." << std::endl;

        return FIELD_POINT_OUTSIDE_DOMAIN;
    }
    if (!b->inside(composite_point)){
        std::cout << "CompositeCurve::composite_field(): composite point " << composite_point << " is outside the boundary." << std::endl;

        return FIELD_POINT_OUTSIDE_DOMAIN;
    }

    // First compute the rarefaction part of the field: dU^-/dxi = r(U^-) for the given family.
    // Notice that the flux and accumulation used are those associated with the
    // rarefaction.
    //
    RealVector reference_vector = composite_object->reference_vector;

    JetMatrix Fm_jet(n);
    rf->jet(rarefaction_point, Fm_jet, 2);
    DoubleMatrix FmJac = Fm_jet.Jacobian();

    JetMatrix Gm_jet(n);
    rg->jet(rarefaction_point, Gm_jet, 2);
    DoubleMatrix GmJac = Gm_jet.Jacobian();

    std::vector<eigenpair> e;
    
    // TODO: All this info may be stored in the rarefaction, since it had to be computed there anyway.
    Eigen::eig(n, FmJac.data(), GmJac.data(), e);

    if (family < 0 || family > e.size() - 1) return FIELD_ERROR;

    // TODO: Check also that lambda is not complex!

    double lambdam = e[family].r;
    RealVector rm(n, e[family].vrr.data()); // Right eigenvector
    RealVector lm(n, e[family].vlr.data()); // Left  eigenvector

    // Verify that rm points in the right direction.
    // Since rm will be used below, it is necessary that this change is permanent.
    //
    if (rm*reference_vector < 0.0) rm = -rm;

    RealVector Gm = Gm_jet.function();

    // For lambda's directional derivative.
    //
    std::vector<DoubleMatrix> Fm_H = Fm_jet.Hessian();
    std::vector<DoubleMatrix> Gm_H = Gm_jet.Hessian();

    RealVector h(Fm_H.size());

    for (int i = 0; i < Fm_H.size(); i++) h(i) = rm*((Fm_H[i] - lambdam*Gm_H[i])*rm);
    
    double dirdrv = (lm*h)/(lm*(GmJac*rm));

    // Now compute the composite part of the field.
    // Notice that the flux and accumulation used are those associated with the
    // composite.
    //
    JetMatrix Fp_jet(n);
    f->jet(composite_point, Fp_jet, 1);
    
    JetMatrix Gp_jet(n);
    g->jet(composite_point, Gp_jet, 1);
    RealVector Gp = Gp_jet.function();

    DoubleMatrix characteristic_matrix = Fp_jet.Jacobian() - lambdam*Gp_jet.Jacobian();

//    if (std::abs(det(characteristic_matrix)) < composite_object->tolerance){
//        std::cout << "Composite field. det = " << det(characteristic_matrix) << std::endl;

//        return FIELD_ERROR; // TODO: Check this tolerance!
//    }

    //
    RealVector cf;
//    int info_solve = solve(characteristic_matrix, dirdrv*(Gp - Gm), cf);

    int info_solve;

    double det_char_matrix = characteristic_matrix(0, 0)*characteristic_matrix(1, 1) - characteristic_matrix(0, 1)*characteristic_matrix(1, 0);

    if (composite_object->compute_first_determinant){
        composite_object->compute_first_determinant = false;
        composite_object->first_determinant = det_char_matrix;
    }

    DoubleMatrix adjugate(2, 2);
    adjugate(0, 0) =  characteristic_matrix(1, 1);
    adjugate(0, 1) = -characteristic_matrix(0, 1);
    adjugate(1, 0) = -characteristic_matrix(1, 0);
    adjugate(1, 1) =  characteristic_matrix(0, 0);

//    cf = adjugate*dirdrv*(Gp - Gm)/det_char_matrix;
    cf = adjugate*dirdrv*(Gp - Gm);


//    if (info_solve == REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR){
//        std::cout << "Composite flux. info_solve = " << info_solve;
//        return FIELD_ERROR;
//    }

    // The first part of the field:
    //
    for (int i = 0; i < n; i++) field[i] = rm(i)*det_char_matrix/composite_object->first_determinant;
//    for (int i = 0; i < n; i++) field[i] = rm(i);

    // The last part of the field:
     //
    for (int i = 0; i < n; i++) field[i + n] = cf(i)/composite_object->first_determinant;
//    for (int i = 0; i < n; i++) field[i + n] = cf(i)/det_char_matrix;

    double nrm = 0.0;
    for (int i = 0; i < n; i++) nrm += cf(i)*cf(i);
    nrm = sqrt(nrm);

    // TODO: This normalization does not appear to improve the quality of the solutions.
    //       The curves are longer, though, which produces better-looking curves.
    //
    //for (int i = 0; i < (*two_n); i++) field[i]first_backwards_rarefaction_point /= nrm;

    // TODO: This normalization seems to be working ok in terms of results. Nevertheless, since
    //       the distance between two consecutive points in the rarefaction thus computed is smaller
    //       than in the original curve, we will not reach the starting point (because).
    //        
    // std::cout << "*** Composite. Field = " << RealVector(*two_n, field) << std::endl;

    // TODO: The last part of the field, corresponding to the composite, MAY need to be normalized.
    //       In that case the first part should be renormalized.

    return FIELD_OK;
}

// TODO: Check with Helmut.
int CompositeCurve::composite_field_near_double_contact3D2D(int *two_n, double *xi, double *pointpair, double *field, int *obj, double* /* Not used */){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    // Extract some information about the composite proper:
    //
    int family = composite_object->family;

    const FluxFunction         *f = composite_object->flux;
    const AccumulationFunction *g = composite_object->accum;
    const Boundary             *b = composite_object->boundary;

    // Extract some information about the rarefaction:
    //
    const FluxFunction         *rf = composite_object->rarflux;
    const AccumulationFunction *rg = composite_object->raraccum;
    const Boundary             *rb = composite_object->rarboundary;

    int n = (*two_n)/2;

    // Extract the points and check if they are within the boundary.
    //
    RealVector rarefaction_point(n);
    RealVector composite_point(n);

    for (int i = 0; i < n; i++){
        rarefaction_point(i) = pointpair[i];
        composite_point(i)   = pointpair[i + n];
    }

    if (!rb->inside(rarefaction_point)){
        std::cout << "CompositeCurve::composite_field(): rarefaction point " << rarefaction_point << " is outside the boundary." << std::endl;

        return FIELD_POINT_OUTSIDE_DOMAIN;
    }
    if (!b->inside(composite_point)){
        std::cout << "CompositeCurve::composite_field(): composite point " << composite_point << " is outside the boundary." << std::endl;

        return FIELD_POINT_OUTSIDE_DOMAIN;
    }

    // First compute the rarefaction part of the field: dU^-/dxi = r(U^-) for the given family.
    // Notice that the flux and accumulation used are those associated with the
    // rarefaction.
    //
    RealVector reference_vector = composite_object->reference_vector;

    JetMatrix Fm_jet(n);
    rf->jet(rarefaction_point, Fm_jet, 2);
    DoubleMatrix FmJac = Fm_jet.Jacobian();

    JetMatrix Gm_jet(n);
    rg->jet(rarefaction_point, Gm_jet, 2);
    DoubleMatrix GmJac = Gm_jet.Jacobian();

    std::vector<eigenpair> e;
    
    // TODO: All this info may be stored in the rarefaction, since it had to be computed there anyway.
    Eigen::eig(n, FmJac.data(), GmJac.data(), e);

    if (family < 0 || family > e.size() - 1) return FIELD_ERROR;

    // TODO: Check also that lambda is not complex!

    double lambdam = e[family].r;
    RealVector rm(n, e[family].vrr.data()); // Right eigenvector
    RealVector lm(n, e[family].vlr.data()); // Left  eigenvector

    // Verify that rm points in the right direction.
    // Since rm will be used below, it is necessary that this change is permanent.
    //
    if (rm*reference_vector < 0.0) rm = -rm;

    // For lambda's directional derivative.
    //
    std::vector<DoubleMatrix> Fm_H = Fm_jet.Hessian();
    std::vector<DoubleMatrix> Gm_H = Gm_jet.Hessian();

    RealVector h(Fm_H.size());

    for (int i = 0; i < Fm_H.size(); i++) h(i) = rm*((Fm_H[i] - lambdam*Gm_H[i])*rm);
    
    double dirdrv = (lm*h)/(lm*(GmJac*rm));

    // Now compute the composite part of the field.
    // Notice that the flux and accumulation used are those associated with the
    // composite.
    //
    JetMatrix Fp_jet(n);
    f->jet(composite_point, Fp_jet, 1);
    
    JetMatrix Gp_jet(n);
    g->jet(composite_point, Gp_jet, 1);

    // Check this part with Helmut.

    RealVector Fp = Fp_jet.function();
    RealVector Gp = Gp_jet.function();

    RealVector Fm = Fm_jet.function();
    RealVector Gm = Gm_jet.function();

    RealVector chi_p(3);
    chi_p(0) = Fp(0)*(Gp(1) - Gm(1)) - Fp(1)*(Gp(0) - Gm(0));
    chi_p(1) = Fp(0)*(Gp(2) - Gm(2)) - Fp(2)*(Gp(0) - Gm(0));
    chi_p(2) = Fp(1)*(Gp(2) - Gm(2)) - Fp(2)*(Gp(1) - Gm(1));

    RealVector chi_m(3);
    chi_m(0) = Fm(0)*(Gp(1) - Gm(1)) - Fm(1)*(Gp(0) - Gm(0));
    chi_m(1) = Fm(0)*(Gp(2) - Gm(2)) - Fm(2)*(Gp(0) - Gm(0));
    chi_m(2) = Fm(1)*(Gp(2) - Gm(2)) - Fm(2)*(Gp(1) - Gm(1));

    double A_correction = (chi_p*chi_m)/(chi_p*chi_p);

    // Characteristic matrix.
    //
    DoubleMatrix cm = A_correction*Fp_jet.Jacobian() - lambdam*Gp_jet.Jacobian();

    //
    RealVector cf;

    int info_solve;

    double det_char_matrix = det(cm);

    if (composite_object->compute_first_determinant){
        composite_object->compute_first_determinant = false;
        composite_object->first_determinant = det_char_matrix;
    }

    DoubleMatrix cofactor(3, 3);

    cofactor(0, 0) =  (cm(1, 1)*cm(2, 2) - cm(1, 2)*cm(2, 1));
    cofactor(0, 1) = -(cm(1, 0)*cm(2, 2) - cm(2, 0)*cm(1, 2));
    cofactor(0, 2) =  (cm(1, 0)*cm(2, 1) - cm(2, 0)*cm(1, 1));

    cofactor(1, 0) = -(cm(0, 1)*cm(2, 2) - cm(2, 1)*cm(0, 2));
    cofactor(1, 1) =  (cm(0, 0)*cm(2, 2) - cm(2, 0)*cm(0, 2));
    cofactor(1, 2) = -(cm(0, 0)*cm(2, 1) - cm(2, 0)*cm(0, 1));

    cofactor(2, 0) =  (cm(0, 1)*cm(1, 2) - cm(1, 1)*cm(0, 2));
    cofactor(2, 1) = -(cm(0, 0)*cm(1, 2) - cm(1, 0)*cm(0, 2));
    cofactor(2, 2) =  (cm(0, 0)*cm(1, 1) - cm(1, 0)*cm(0, 1));


    DoubleMatrix adjugate(transpose(cofactor));

    double det_cm = det(cm);

//    cf = adjugate*dirdrv*(Gp - Gm)/det_char_matrix;
    cf = adjugate*dirdrv*(Gp - Gm);


//    if (info_solve == REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR){
//        std::cout << "Composite flux. info_solve = " << info_solve;
//        return FIELD_ERROR;
//    }

    // The first part of the field:
    //
    for (int i = 0; i < n; i++) field[i] = rm(i)*det_char_matrix/composite_object->first_determinant;
//    for (int i = 0; i < n; i++) field[i] = rm(i);

    // The last part of the field:
     //
    for (int i = 0; i < n; i++) field[i + n] = cf(i)/composite_object->first_determinant;
//    for (int i = 0; i < n; i++) field[i + n] = cf(i)/det_char_matrix;

    double nrm = 0.0;
    for (int i = 0; i < n; i++) nrm += cf(i)*cf(i);
    nrm = sqrt(nrm);

    // TODO: This normalization does not appear to improve the quality of the solutions.
    //       The curves are longer, though, which produces better-looking curves.
    //
    //for (int i = 0; i < (*two_n); i++) field[i]first_backwards_rarefaction_point /= nrm;

    // TODO: This normalization seems to be working ok in terms of results. Nevertheless, since
    //       the distance between two consecutive points in the rarefaction thus computed is smaller
    //       than in the original curve, we will not reach the starting point (because).
    //        
    // std::cout << "*** Composite. Field = " << RealVector(*two_n, field) << std::endl;

    // TODO: The last part of the field, corresponding to the composite, MAY need to be normalized.
    //       In that case the first part should be renormalized.

    return FIELD_OK;
}

RealVector CompositeCurve::composite_field(const RealVector &final_point_pair){

    double xi = 0.0;
    
    int two_n = final_point_pair.size();
    int n     = two_n/2;
    
    RealVector field(two_n);
    
    int info = composite_field(&two_n, &xi, (RealVector(final_point_pair)).components(), field.components(), (int*)this, 0);

//    std::cout << "CompositeCurve. Field at the last point = " << field << std::endl;

    return RealVector(n, &(field.components()[n]));
}

void CompositeCurve::all_eigenvalues(const RealVector &p, int family, RealVector &point_eigenvalues){
    std::vector<double> lambda;
    Eigen::fill_eigenvalues(flux, accum, p, lambda);

    point_eigenvalues.resize(lambda.size());
    for (int i = 0; i < lambda.size(); i++) point_eigenvalues(i) = lambda[i];

    return;
}

void CompositeCurve::add_point_to_curve(const RealVector &p, int back, const Curve &rarcurve, Curve &curve){
    curve.curve.push_back(p);
    curve.back_pointer.push_back(back);

    RealVector point_eigenvalues;
    all_eigenvalues(p, family, point_eigenvalues);

    curve.eigenvalues.push_back(point_eigenvalues);

//    curve.speed.push_back(point_eigenvalues(family));
    curve.speed.push_back(rarcurve.speed[back]);

    

    return;
}

int CompositeCurve::curve(const AccumulationFunction *RarAccum, const FluxFunction *RarFlux,
                          const Boundary *RarBoundary, 
                          const Curve &rarcurve,
                          const RealVector &composite_initial_point,
                          int last_point_in_rarefaction,
                          const ODE_Solver *odesolver,
                          double deltaxi,
                          void *linobj, double (*linear_function)(void *o, const RealVector &p),
                          int where_composite_begins, int fam, 
                          Curve &new_rarcurve,
                          Curve &compositecurve,
                          RealVector &final_direction,
                          int &reason_why,
                          int &edge){

    normalize_with_respect_to_whom = NORMALIZE_WITH_RESPECT_TO_RAREFACTION;
    compute_first_determinant = true;
    cmp_deltaxi = deltaxi;
    use_field_near_double_contact = false;

    int n     = rarcurve.curve[0].size();
    int two_n = 2*n;

    family = fam;
    lambda_init_base_rarefaction = rarcurve.speed[0]; // TODO: Maybe this will not be the zeroeth element.

    compositecurve.clear();
    compositecurve.type = COMPOSITE_CURVE;
    compositecurve.family = family;

    new_rarcurve.clear();
    new_rarcurve.type = RAREFACTION_CURVE;
    new_rarcurve.family = family;

    // To be used by the field:
    //
    rarflux     = RarFlux;
    raraccum    = RarAccum;
    rarboundary = RarBoundary;

    RealVector rarcmp_point(two_n);

    // Store the first point.
    //
    add_point_to_curve(composite_initial_point, last_point_in_rarefaction, rarcurve, compositecurve);
    new_rarcurve.curve.push_back(rarcurve.curve[last_point_in_rarefaction]);

    // Initialize the intersection of the rarefaction and a line.
    //
    double old_linear_function_value = 0.0;
    if (linear_function != 0) old_linear_function_value = (*linear_function)(linobj, composite_initial_point);
    double linear_function_value = old_linear_function_value;

    // Initialize the composite curve.
    //
    int index_of_corresponding_point_in_rarefaction;

    if (where_composite_begins == COMPOSITE_BEGINS_AT_INFLECTION){
        std::cout << "CompositeCurve: COMPOSITE_BEGINS_AT_INFLECTION" << std::endl;

        // TODO: This retreat may be insufficient for some cases. What to do: if the auxiliary shockcurve reaches the boundary,
        //       which is evidently wrong, increase the value of retreat and try again.
        int retreat = 3;
    
        int index_last_rar = rarcurve.curve.size() - 1; // TODO: This may be replaced by an index, not necessarily the last element of the rarefaction will be used.
        RealVector rar_last_point = rarcurve.curve[index_last_rar]; // TODO: See line above.

        // Retreat some points from the rarefaction's end at inflection.
        if (rarcurve.curve.size() > retreat){
            // Find the point beyond the inflection curve
            int retreat_index = index_last_rar - retreat; // rarcurve.size() - retreat - 1;
            RealVector first_backwards_rarefaction_point = rarcurve.curve[retreat_index];
            
            RealVector direction = rar_last_point - first_backwards_rarefaction_point;
            RealVector p = rar_last_point + 2.0*direction; // This is the initial point for the Newton iteration which corresponds to the retreat point of

            // Compute the Hugoniot locus that with first_backwards_rarefaction_point as reference point
            // Use the Newton method to find a point on the Composite curve corresponding to the rarefaction going back by (retreat). 
            //
            ReferencePoint ref(first_backwards_rarefaction_point, RarFlux, RarAccum, (Viscosity_Matrix*)0);
            shock->set_reference_point(ref);

            RealVector first_composite_point;
            int info = shock->find_point_for_sigma_equal_reference_lambda(p, rarcurve.speed[retreat_index], first_composite_point);

            if (info == SHOCKCURVE_NEWTON_DID_NOT_CONVERGE){
                // In this case, try to retreat even more and create a new shockcurve.
                Curve shockcurve;
                std::vector<int> transition_current_index, transition_current_family, transition_reference_index, transition_reference_family;
                int edge;

                //retreat = 20;

                RealVector initial_point_for_shock(rarcurve.curve[rarcurve.curve.size() - retreat]);

                int shock_stopped_because; // TODO: Will this be used???
                int info_shock = shock->curve_engine(ReferencePoint(initial_point_for_shock, RarFlux, RarAccum, 0), initial_point_for_shock, 
                                                     initial_point_for_shock - rarcurve.curve[rarcurve.curve.size() - retreat - 1], family, 
                                                     SHOCKCURVE_SHOCK_CURVE, 
                                                     SHOCK_SIGMA_EQUALS_LAMBDA_OF_FAMILY_AT_LEFT, DONT_CHECK_EQUALITY_AT_RIGHT,
                                                     USE_CURRENT_FAMILY,
                                                     STOP_AFTER_TRANSITION,
                                                     shockcurve, 
                                                     transition_current_index,
                                                     transition_current_family,
                                                     transition_reference_index,
                                                     transition_reference_family,  
                                                     shock_stopped_because,
                                                     edge);

                if (info_shock == SHOCKCURVE_OK){
                    first_backwards_rarefaction_point = initial_point_for_shock;
                    first_composite_point = shockcurve.curve.back();
                }
                else return COMPOSITE_ERROR;
            }

            if (!boundary->inside(first_composite_point)){
                // TODO: Maybe boundary->intersection() should be used. Decide it.
                std::cout << "Out was outside the domain!" << std::endl;
                return COMPOSITE_ERROR_AT_BEGINNING_OUT_OF_BOUNDARY;
            }

            // TODO: If the composite does not start at the inflection, the last valid index of the rarefaction must be passed.
            //
            reference_vector = first_backwards_rarefaction_point - rarcurve.curve.back();
            normalize(reference_vector);

            for (int i = 0; i < n; i++){
                rarcmp_point(i)     = first_backwards_rarefaction_point(i);
                rarcmp_point(i + n) = first_composite_point(i);
            }

            add_point_to_curve(RealVector(n, n, rarcmp_point), retreat_index, rarcurve, compositecurve);

            // The second point in the composite corresponds to this point in the rarefaction:
            //
            index_of_corresponding_point_in_rarefaction = index_last_rar - retreat;
        }
        else { // rarcurve.curve.size() > retreat
        }
    }
    else if (where_composite_begins == COMPOSITE_AFTER_COMPOSITE){
        std::cout << "CompositeCurve: COMPOSITE_AFTER_COMPOSITE" << std::endl;

        // TODO: If the composite does not start at the inflection, the last valid index of the rarefaction must be passed.
        //
        reference_vector = rarcurve.curve[rarcurve.curve.size() - 2] - rarcurve.curve[rarcurve.curve.size() - 1];
        normalize(reference_vector);

        for (int i = 0; i < n; i++){
            rarcmp_point(i)     = rarcurve.curve[last_point_in_rarefaction /*rarcurve.curve.size() - 1*/](i); //first_backwards_rarefaction_point(i);
            rarcmp_point(i + n) = composite_initial_point(i);
        }

        // The second point in the composite corresponds to this point in the rarefaction:
        //
        index_of_corresponding_point_in_rarefaction = last_point_in_rarefaction - 1;
    }

    // Compute the rest of the composite curve.

    double init_time = 0.0;
    double final_time = init_time + deltaxi;

    // Out
    RealVector out;

    // The sign of the determinant of the characteristic matrix changes when the composite passes through the double contact.
    // To detect this event, use the variables below.
    //
    // Compute the determinant before proceeding (at the first point of the rarefaction-composite subspace).
    //
    double current_determinant, previous_determinant;
    int info_previous_determinant = double_contact_signal_event(rarcmp_point, previous_determinant, (int*)this, 0);

    if (info_previous_determinant == BISECTION_FUNCTION_ERROR){
        std::cout << "Composite. Error = " << COMPOSITE_ERROR_AT_BEGINNING_DETERMINANT << ". Error in file \""<< __FILE__ << "\", method integrated() at line " << __LINE__ << std::endl;
            
        return COMPOSITE_ERROR_AT_BEGINNING_DETERMINANT;
    }

    // The sign of (sigma - lambda_ref) changes when the composite consumes the rarefaction from whence it spawned.
    // To detect this event, use the variables below.
    //
    // Compute (sigma - lambda_ref).
    //
    double previous_diff_lambda_init, current_diff_lambda_init;
    int info_previous_diff_lambda_init = rarefaction_of_composite_signal_event(rarcmp_point, previous_diff_lambda_init, (int*)this, 0);

    field = &composite_field;

    while (true){
//        std::cout << "CompositeCurve, inside while. compositecurve.curve.size() = " << compositecurve.curve.size() << std::endl;

        int info_odesolver = odesolver->integrate_step(field, (int*)this, (double*)0 /*function_data*/, 
                                                       init_time,  rarcmp_point,
                                                       final_time, out);

        if (use_field_near_double_contact){
            field = &composite_field_near_double_contact;
        }

//        int info_odesolver = odesolver->integrate_step(&composite_field, (int*)this, (double*)0 /*function_data*/, 
//                                                       init_time,  rarcmp_point,
//                                                       final_time, out);

//        int info_odesolver = odesolver->integrate_step(&composite_field_near_double_contact, (int*)this, (double*)0 /*function_data*/, 
//                                                       init_time,  rarcmp_point,
//                                                       final_time, out);

        // Update the index of the corresponding point in the rarefaction.
        //
        index_of_corresponding_point_in_rarefaction--;


        if (info_odesolver == ODE_SOLVER_ERROR){
            std::cout << "CompositeCurve. info_odesolver == ODE_SOLVER_ERROR." << std::endl;
            return COMPOSITE_ERROR;
        }
                
        // Has the composite curve reached the boundary?
        //
        RealVector point_on_composite(n, n, out);
        RealVector prev_point_on_composite(n, n, rarcmp_point);
        RealVector intersection_point;
        int info_intersect = boundary->intersection(point_on_composite, prev_point_on_composite, intersection_point, edge);

        // TODO: Remove from here and place it immediately after the call to odesolver.
        // What to do with the corresponding point in the rarefaction curve?
        // Right now the composite ends in a point that is unmatched in the rarefaction.
        //
        if (info_intersect == BOUNDARY_INTERSECTION_FOUND){
            compositecurve.last_point = intersection_point;
            add_point_to_curve(intersection_point, index_of_corresponding_point_in_rarefaction, rarcurve, compositecurve);
            new_rarcurve.curve.push_back(RealVector(0, n, rarcmp_point)); // TODO: WRONG!!!

            reason_why = COMPOSITE_REACHED_BOUNDARY;
            compositecurve.reason_to_stop = COMPOSITE_REACHED_BOUNDARY;

            std::cout << "Composite reached the boundary" << std::endl;

            return COMPOSITE_OK;
        }

        // Stop criteria, for the time being: if det(characteristic_matrix) changes sign, stop (and perform a linear interpolation).
        // As a more refined strategy, create a signal_event watching (lamba_p - lambda_m).

        // Signal event
        int info_current_determinant = double_contact_signal_event(out, current_determinant, (int*)this, 0);

        if (info_current_determinant == BISECTION_FUNCTION_ERROR){
            std::cout << "Composite. Error = " << COMPOSITE_ERROR_AT_DETERMINANT << ". Error in file \"" << __FILE__ << "\", method integrated() at line " << __LINE__ << std::endl;
            
            compositecurve.reason_to_stop = COMPOSITE_ERROR_AT_DETERMINANT;
            return COMPOSITE_ERROR_AT_DETERMINANT;
        }
                
        if (current_determinant*previous_determinant <= 0.0){
            std::cout << "Composite. Will invoke Bisection (Determinant changed sign)." << std::endl;

            double bisection_epsilon = 1e-20;
                    
            // Output here:
            double c_t;
            RealVector p_c;
                
//            int info_bisection = Bisection::bisection_method(init_time,  rarcmp_point,
//                                                             final_time, out,
//                                                             bisection_epsilon, 
//                                                             c_t, p_c,
//                                                             &composite_field, (int*)this, (double*)0,
//                                                             odesolver, 
//                                                             &double_contact_signal_event, (int*)this /*int *signal_event_object*/, 0 /*int *signal_event_data*/);

            int info_bisection = Bisection::bisection_method(init_time,  rarcmp_point,
                                                             final_time, out,
                                                             bisection_epsilon, 
                                                             c_t, p_c,
                                                             field, (int*)this, (double*)0,
                                                             odesolver, 
                                                             &double_contact_signal_event, (int*)this /*int *signal_event_object*/, 0 /*int *signal_event_data*/);

            // Return the last point and the final direction.
            if (info_bisection == BISECTION_FUNCTION_OK){
                compositecurve.last_point = RealVector(n, n, p_c);
                add_point_to_curve(compositecurve.last_point, index_of_corresponding_point_in_rarefaction, rarcurve, compositecurve);
                new_rarcurve.curve.push_back(RealVector(0, n, p_c));

                compositecurve.final_direction = final_direction = RealVector(n, n, out - rarcmp_point); // composite_field(p_c);
                normalize(final_direction);
                normalize(compositecurve.final_direction);
                        
                reason_why = COMPOSITE_REACHED_DOUBLE_CONTACT;
                compositecurve.reason_to_stop = COMPOSITE_REACHED_DOUBLE_CONTACT;

                std::cout << "Composite will end now (double contact). final_direction = " << final_direction << std::endl;

                return COMPOSITE_OK;
            }
            else {
                std::cout << "Composite. Error = " << COMPOSITE_ERROR << ". Error in file \""<< __FILE__ << "\", method integrated() at line " << __LINE__ << std::endl;

                compositecurve.reason_to_stop = COMPOSITE_ERROR_AT_DETERMINANT;
                return COMPOSITE_ERROR_AT_DETERMINANT;
            }
        }

        if (linear_function != 0){
            linear_function_value = (*linear_function)(linobj, point_on_composite);

            if (linear_function_value*old_linear_function_value < 0.0){
                double alpha = linear_function_value/(linear_function_value - old_linear_function_value);

                RealVector point_on_line = alpha*prev_point_on_composite + (1.0 - alpha)*point_on_composite;

                add_point_to_curve(point_on_line, index_of_corresponding_point_in_rarefaction, rarcurve, compositecurve);

                reason_why = COMPOSITE_REACHED_LINE;
                compositecurve.reason_to_stop = COMPOSITE_REACHED_LINE;

                compositecurve.last_point = point_on_line;
                compositecurve.final_direction = point_on_composite - prev_point_on_composite; //rarcurve.curve.back() - rarcurve.curve[rarcurve.curve.size() - 2];
                normalize(compositecurve.final_direction);

                final_direction = compositecurve.final_direction;

                return COMPOSITE_OK;
            }
            else old_linear_function_value = linear_function_value;
        }
                
        // Check if the beginning of the rarefaction was reached.
        int info_current_diff_lambda_init = rarefaction_of_composite_signal_event(out, current_diff_lambda_init, (int*)this, 0);
        if (info_current_diff_lambda_init == BISECTION_FUNCTION_ERROR){

            compositecurve.reason_to_stop = COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING;
            return COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING;
        }
            

//        std::cout << "CompositeCurve: previous_diff_lambda_init = " << previous_diff_lambda_init << ", current_diff_lambda_init = " << current_diff_lambda_init << std::endl;
        //TestTools::pause();
   
        if (previous_diff_lambda_init*current_diff_lambda_init <= 0.0){
            std::cout << "Composite. Near the beginning of the rarefaction." << std::endl;

            // Bisection here.
            double bisection_epsilon = 1e-8;
                    
            double c_t;
            RealVector p_c;
                    
            int info_bisection = Bisection::bisection_method(init_time,  rarcmp_point,
                                                             final_time, out,
                                                             bisection_epsilon, 
                                                             c_t, p_c,
                                                             &composite_field, (int*)this, (double*)0,
                                                             odesolver,
                                                             &rarefaction_of_composite_signal_event, (int*)this /*int *signal_event_object*/, 0 /*int *signal_event_data*/);
                                                                     
            if (info_bisection == BISECTION_FUNCTION_OK){
                std::cout << "out = " << out << ", p_c = " << p_c << std::endl;

                compositecurve.last_point = RealVector(n, n, p_c);
                add_point_to_curve(compositecurve.last_point, index_of_corresponding_point_in_rarefaction, rarcurve, compositecurve);
                new_rarcurve.curve.push_back(RealVector(0, n, p_c));
                    
//                compositecurve.final_direction = final_direction = composite_field(p_c);
                compositecurve.final_direction = final_direction = RealVector(n, n, out - rarcmp_point); // composite_field(p_c);
                normalize(final_direction);
                normalize(compositecurve.final_direction);
                    
                reason_why = COMPOSITE_COMPLETED;
                compositecurve.reason_to_stop = COMPOSITE_COMPLETED;

                std::cout << "Composite will end now (composite completed). final_direction = " << final_direction << std::endl;
                return COMPOSITE_OK;
            }
            else {
                std::cout << "Composite. Near the beginning of the rarefaction. Bisection error!" << std::endl;

                compositecurve.reason_to_stop = COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING;
                return COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING;
            }
        }

        // Explicit bifurcation curve intersection. TEMPORARY!
        if (explicit_bifurcation_curve != 0){
            normalize_with_respect_to_whom = NORMALIZE_WITH_RESPECT_TO_COMPOSITE;

            int info_transition = transition_with_explicit_bifurcation(odesolver, rarcmp_point, init_time, out, final_time);

            if (info_transition == SECUNDARY_BIFURCATION_DETECTED){
                std::cout << "CompositeCurve. Transition detected @ " << out << std::endl;
                compositecurve.explicit_bifurcation_transition_index.push_back(compositecurve.curve.size());
            }

            normalize_with_respect_to_whom = NORMALIZE_WITH_RESPECT_TO_RAREFACTION;
        }
        // Explicit bifurcation curve intersection. TEMPORARY!
                
        // Update the determinant.
        //
        previous_determinant = current_determinant;

        //Update the lambda_diff.
        //
        previous_diff_lambda_init = current_diff_lambda_init;

        add_point_to_curve(point_on_composite, index_of_corresponding_point_in_rarefaction, rarcurve, compositecurve);
        new_rarcurve.curve.push_back(RealVector(0, n, rarcmp_point));

        // Update.
        //
        for (int i = 0; i < n; i++) reference_vector(i) = out(i) - rarcmp_point(i);
        normalize(reference_vector);

        init_time = final_time;
        final_time += deltaxi;

        rarcmp_point = out;
    } // while
}


// This method should be used after CompositeCurve::curve() is called and what comes after it is a COMPOSITE_AFTER_COMPOSITE.
// It is assumed that wavecurve comes with enough information.
//
int CompositeCurve::correct_last_point(const ODE_Solver *odesolver, double deltaxi, WaveCurve &wavecurve){
    // First, the value of lambda must be found.

    // Find the rarefaction that came before the last one.
    // The last curve is a composite, the next-to-last is the rarefaction that was exhausted.
    // Skip the last two curves and find the next rarefaction in reverse order.
    //
    int rarefaction_to_be_used = wavecurve.wavecurve.size() - 2; 

    while (true){
        rarefaction_to_be_used--;

        if (rarefaction_to_be_used < 0){
            std::cout << "CompositeCurve::correct_last_point(): The wavecurve is too short! Aborting." << std::endl;

            return COMPOSITE_LAST_POINT_ERROR;
        }
        else if (wavecurve.wavecurve[rarefaction_to_be_used].type == RAREFACTION_CURVE) break;
    }

    // The composite curve that comes after the abovementioned rarefaction (because the last point in the composite 
    // is related to the point in the rarefaction whose lambda is needed).
    //
    int composite_to_be_used = rarefaction_to_be_used + 1;
    int last_point_in_rarefaction_to_be_used = wavecurve.wavecurve[composite_to_be_used].back_pointer.back();

    // This is the last used point of the rarefaction.
    //
    RealVector extremum_rarefaction_point = wavecurve.wavecurve[rarefaction_to_be_used].curve[last_point_in_rarefaction_to_be_used];

    // The dimension of the space.
    //
    int n = extremum_rarefaction_point.size(), two_n = 2*n;

    // The flux and accumulation at this extremum will be used (potentially) several times to compute
    // the speed. Store them.
    //
    JetMatrix rar_F_jet(n), rar_G_jet(n);
    rarflux->jet(extremum_rarefaction_point, rar_F_jet, 0);
    raraccum->jet(extremum_rarefaction_point, rar_G_jet, 0);

    RealVector rarF = rar_F_jet.function();
    RealVector rarG = rar_G_jet.function();

    // This is the value of lambda at the last used point of the rarefaction.
    //
    double lambda = wavecurve.wavecurve[rarefaction_to_be_used].speed[last_point_in_rarefaction_to_be_used];

    std::cout << "CompositeCurve::correct_last_point(). last_point_in_previous_rarefaction = " << last_point_in_rarefaction_to_be_used << std::endl;

    // Second, the rarefaction-composite field must be integrated, until a point is found where the sign of (sigma - lambda) changes.
    // Here: lambda is the abovementioned value. Sigma is the speed between the point in the rarefaction curve whose speed is lambda and
    // the current point in the composite curve.
    //
    // The integration will commence from the last point of the wavecurve (in its current state). The value
    // of (sigma - lambda) will be initialized accordingly.

    // The point in the rarefaction-composite space.
    //
    int corresponding_index_in_rarefaction = wavecurve.wavecurve.back().back_pointer.back();
    int corresponding_rarefaction = wavecurve.wavecurve.back().back_curve_index;

    RealVector previous_composite_point   = wavecurve.wavecurve.back().curve.back();
    RealVector previous_rarefaction_point = wavecurve.wavecurve[corresponding_rarefaction].curve[corresponding_index_in_rarefaction];

    RealVector previous_point(two_n);

    for (int i = 0; i < n; i++){
        previous_point(i)     = previous_rarefaction_point(i);
        previous_point(i + n) = previous_composite_point(i);
    } 

    std::cout << "Combined point: " << previous_point << std::endl;
    std::cout << "    Final direction: " << wavecurve.wavecurve.back().final_direction << std::endl;

    // The output of the iteration.
    //
    RealVector point(previous_point);

    double init_time = 0.0;

    // The sign of (sigma - lambda) should change at some point (hopefully, after only one iteration...)
    //
    RealVector composite_point(n, n, point);

    JetMatrix cmp_F_jet(n), cmp_G_jet(n);
    flux->jet(composite_point, cmp_F_jet, 0);
    accum->jet(composite_point, cmp_G_jet, 0);

    double previous_sigma_minus_lambda = shock->get_HugoniotContinuation()->sigma(rarF, rarG, cmp_F_jet.function(), cmp_G_jet.function()) - lambda;

    reference_vector = wavecurve.wavecurve[corresponding_rarefaction].curve[corresponding_index_in_rarefaction + 1] - previous_rarefaction_point;
    normalize(reference_vector);

    int it = 0;

    while (true){
        it++;
        std::cout << "Correction. it = " << it << std::endl;

        // Compute the new point.
        //
        RealVector point;

        double final_time = init_time + deltaxi;

        int info_odesolver = odesolver->integrate_step(&composite_field, (int*)this, (double*)0 /*function_data*/, 
                                                       init_time,  previous_point,
                                                       final_time, point);    

        // Compute the value of (sigma - lambda) in the new point and compare it with the
        // one from the previous iteration.
        //
        RealVector composite_point(n, n, point);

        JetMatrix cmp_F_jet(n), cmp_G_jet(n);
        flux->jet(composite_point, cmp_F_jet, 0);
        accum->jet(composite_point, cmp_G_jet, 0);

        double sigma_minus_lambda = shock->get_HugoniotContinuation()->sigma(rarF, rarG, cmp_F_jet.function(), cmp_G_jet.function()) - lambda;

        if (sigma_minus_lambda*previous_sigma_minus_lambda <= 0.0){
            std::cout << "Number of iterations: " << it << std::endl;
            std::cout << "sigma_minus_lambda = " << sigma_minus_lambda << ", previous_sigma_minus_lambda = " << previous_sigma_minus_lambda << std::endl;

            // Proceed to bisection. 
            //
            double bisection_epsilon = 1e-10; // Problems with this threshold.
                    
            // TODO: These names could change to match the name of the signal event method.
            //
            lambda_at_double_contact = lambda;
            rar_F_at_double_contact = rarF;
            rar_G_at_double_contact = rarG;

            double c_t;
            RealVector p_c;
                    
            int info_bisection = Bisection::bisection_method(init_time,  previous_point,
                                                             final_time, point,
                                                             bisection_epsilon, 
                                                             c_t, p_c,
                                                             &composite_field, (int*)this, (double*)0,
                                                             odesolver,
                                                             &sigma_minus_lambda_signal_event, (int*)this /*int *signal_event_object*/, 0 /*int *signal_event_data*/);

            std::cout << "Correction: p_c = " << p_c << std::endl;

            JetMatrix correction_F_jet(n), correction_G_jet(n);
            flux->jet(RealVector(n, n, p_c), correction_F_jet, 0);
            accum->jet(RealVector(n, n, p_c), correction_G_jet, 0);

            std::cout << "************** = > Sigma = " << shock->get_HugoniotContinuation()->sigma(rarF, rarG, correction_F_jet.function(), correction_G_jet.function()) << ", lambda = " << lambda << std::endl;
            std::cout << "    sigma - lambda = " << shock->get_HugoniotContinuation()->sigma(rarF, rarG, correction_F_jet.function(), correction_G_jet.function()) - lambda << std::endl;
            
            if (info_bisection == BISECTION_FUNCTION_OK){
                // TODO: Replace the last point for the computed by the Bisection.
                std::cout << "$$$$ Before replacing: " << wavecurve.wavecurve.back().curve.back() << std::endl;
                std::cout << "$$$$  After replacing: " << RealVector(n, n, p_c) << std::endl; 

                wavecurve.wavecurve.back().curve.back() = RealVector(n, n, p_c);
                //wavecurve.wavecurve.back().speed.back() = lambda; // segfault here.

                //TestTools::pause("Correction OK!");

                return COMPOSITE_OK;
            }
            else {
                return COMPOSITE_LAST_POINT_ERROR;
            }
        }
        else {
            previous_sigma_minus_lambda = sigma_minus_lambda;
            previous_point = point;
        }

    }

//    double previous_sigma_minus_lambda = sigma - lambda;

//    while (true){
//        double sigma = shock->get_HugoniotContinuation()->shockspeed(rarflux, raraccum, extremum_rarefaction_point,
//                                                                     flux,    accum,    composite_point);

//        std::cout << "   ===> sigma = " << sigma << std::endl;


//        int info_odesolver = odesolver->integrate_step(&composite_field, (int*)this, (double*)0 /*function_data*/, 
//                                                       init_time,  rarcmp_point,
//                                                       final_time, out);    

//        if (previous_diff_lambda_init*current_diff_lambda_init <= 0.0){
//            std::cout << "Composite. Near the beginning of the rarefaction." << std::endl;

//            // Bisection here.
//            double bisection_epsilon = 1e-20;
//                    
//            double c_t;
//            RealVector p_c;
//                    
//            int info_bisection = Bisection::bisection_method(init_time,  rarcmp_point,
//                                                             final_time, out,
//                                                             bisection_epsilon, 
//                                                             c_t, p_c,
//                                                             &composite_field, (int*)this, (double*)0,
//                                                             odesolver,
//                                                             &rarefaction_of_composite_signal_event, (int*)this /*int *signal_event_object*/, 0 /*int *signal_event_data*/);
//                                                                     
//            if (info_bisection == BISECTION_FUNCTION_OK){
//                std::cout << "out = " << out << ", p_c = " << p_c << std::endl;

//                compositecurve.last_point = RealVector(n, n, p_c);
//                add_point_to_curve(compositecurve.last_point, index_of_corresponding_point_in_rarefaction, compositecurve);
//                    
////                compositecurve.final_direction = final_direction = composite_field(p_c);
//                compositecurve.final_direction = final_direction = RealVector(n, n, out - rarcmp_point); // composite_field(p_c);
//                normalize(final_direction);
//                    
//                reason_why = COMPOSITE_COMPLETED;
//                compositecurve.reason_to_stop = COMPOSITE_COMPLETED;

//                std::cout << "Composite will end now (composite completed). final_direction = " << final_direction << std::endl;
//                return COMPOSITE_OK;
//            }
//            else {
//                std::cout << "Composite. Near the beginning of the rarefaction. Bisection error!" << std::endl;

//                compositecurve.reason_to_stop = COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING;
//                return COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING;
//            }
//        }

//    }

    return COMPOSITE_OK;
}

CompositeCurve::CompositeCurve(const AccumulationFunction *a, const FluxFunction *f, const Boundary *b, ShockCurve *s, Explicit_Bifurcation_Curves *ebc){
    flux = f;
    accum = a;
    boundary = b;
    
    shock = s;

    explicit_bifurcation_curve = ebc;

    tolerance = 1e-10; // For the determinant of the characteristic matrix, qv composite_field().


//    retreat = 20; // For Helmut, 20.
    retreat = 2; // Number of points to be skipped from the rarefaction's end when the Composite starts on the inflection curve.
}

CompositeCurve::~CompositeCurve(){
}

// TODO: CHECK THIS!!!!
int CompositeCurve::double_contact_signal_event(const RealVector &where, double & determinant, int *obj, int * /*not used*/){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    // Extract some information about the composite proper:
    //
    int family = composite_object->family;

    const FluxFunction         *f = composite_object->flux;
    const AccumulationFunction *g = composite_object->accum;
    const Boundary             *b = composite_object->boundary;

    // Extract some information about the rarefaction:
    //
    const FluxFunction         *rf = composite_object->rarflux;
    const AccumulationFunction *rg = composite_object->raraccum;
    const Boundary             *rb = composite_object->rarboundary;

    int two_n = where.size();
    int n     = two_n/2;

    // Extract the points and check if they are within the boundary.
    //
    RealVector rarefaction_point(0, n, where);
    RealVector composite_point(n, n, where);

    // This may be overkill. Remove at will.
    //
    if (!rb->inside(rarefaction_point)) return BISECTION_FUNCTION_ERROR;
    if (!b->inside(composite_point))    return BISECTION_FUNCTION_ERROR;

    // First compute the rarefaction part of the field: dU^-/dxi = r(U^-) for the given family.
    // Notice that the flux and accumulation used are those associated with the
    // rarefaction.
    //
    //RealVector reference_vector = composite_object->reference_vector;

    JetMatrix Fm_jet(n);
    rf->jet(rarefaction_point, Fm_jet, 1);
    DoubleMatrix FmJac = Fm_jet.Jacobian();

    JetMatrix Gm_jet(n);
    rg->jet(rarefaction_point, Gm_jet, 1);
    DoubleMatrix GmJac = Gm_jet.Jacobian();

    std::vector<eigenpair> e;
    
    Eigen::eig(n, FmJac.data(), GmJac.data(), e);

    //if (family < 0 || family > e.size() - 1) return FIELD_ERROR;

    // TODO: Check also that lambda is not complex!

    double lambda = e[family].r;

    // Now compute the composite part of the field.
    // Notice that the flux and accumulation used are those associated with the
    // composite.
    //
    JetMatrix Fp_jet(n);
    f->jet(composite_point, Fp_jet, 1);

    JetMatrix Gp_jet(n);
    g->jet(composite_point, Gp_jet, 1);

    determinant = det(Fp_jet.Jacobian() - lambda*Gp_jet.Jacobian());
    
    return BISECTION_FUNCTION_OK;    
}

//int CompositeCurve::characteristic_shock_signal_event(const RealVector &where, double &diff_lambda, int *obj, int * /*not used*/){
//    CompositeCurve *composite_object = (CompositeCurve*)obj;

//    // Extract some information about the composite proper:
//    //
//    int family = composite_object->family;

//    const FluxFunction         *f = composite_object->flux;
//    const AccumulationFunction *g = composite_object->accum;
//    const Boundary             *b = composite_object->boundary;

//    // Extract some information about the rarefaction:
//    //
//    const FluxFunction         *rf = composite_object->rarflux;
//    const AccumulationFunction *rg = composite_object->raraccum;
//    const Boundary             *rb = composite_object->rarboundary;

//    int two_n = where.size();
//    int n     = two_n/2;

//    // Extract the points and check if they are within the boundary.
//    //
//    RealVector rarefaction_point(0, n, where);
//    RealVector composite_point(n, n, where);

//    if (!rb->inside(rarefaction_point)) return BISECTION_FUNCTION_ERROR;
//    if (!b->inside(composite_point))    return BISECTION_FUNCTION_ERROR;

//    // Compute lambda-
//    //
//    JetMatrix Fm_jet(n);
//    rf->jet(rarefaction_point, Fm_jet, 1);
//    DoubleMatrix FmJac = Fm_jet.Jacobian();

//    JetMatrix Gm_jet(n);
//    rg->jet(rarefaction_point, Gm_jet, 1);
//    DoubleMatrix GmJac = Gm_jet.Jacobian();

//    std::vector<eigenpair> em;
//    
//    Eigen::eig(n, FmJac.data(), GmJac.data(), em); 
// 
//    double lambdam = em[family].r;  

//    // Compute lambda+
//    //
//    JetMatrix Fp_jet(n);
//    f->jet(composite_point, Fp_jet, 1);
//    DoubleMatrix FpJac = Fp_jet.Jacobian();

//    JetMatrix Gp_jet(n);
//    g->jet(composite_point, Gp_jet, 1);
//    DoubleMatrix GpJac = Gp_jet.Jacobian();

//    std::vector<eigenpair> ep;
//    
//    Eigen::eig(n, FpJac.data(), GpJac.data(), ep); 
// 
//    double lambdap = ep[family].r;  

//    diff_lambda = lambdap - lambdam;

//    return BISECTION_FUNCTION_OK;
//}

int CompositeCurve::rarefaction_of_composite_signal_event(const RealVector &where, double & current_diff_lambda, int *obj, int * /*not used*/){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    // Extract some information about the composite proper:
    //
    int family = composite_object->family;
    double lambda_init_base_rarefaction = composite_object->lambda_init_base_rarefaction;

    // Extract some information about the rarefaction:
    //
    const FluxFunction         *rf = composite_object->rarflux;
    const AccumulationFunction *rg = composite_object->raraccum;
    const Boundary             *rb = composite_object->rarboundary;

    int two_n = where.size();
    int n     = two_n/2;

    // Extract the points and check if they are within the boundary.
    //
    RealVector rarefaction_point(n);
    
    for (int i = 0; i < n; i++) rarefaction_point(i) = where(i);

    if (!rb->inside(rarefaction_point)) return BISECTION_FUNCTION_ERROR;
    
    JetMatrix Fm_jet(n);
    rf->jet(rarefaction_point, Fm_jet, 1);
    DoubleMatrix FmJac = Fm_jet.Jacobian();

    JetMatrix Gm_jet(n);
    rg->jet(rarefaction_point, Gm_jet, 1);
    DoubleMatrix GmJac = Gm_jet.Jacobian();

    std::vector<eigenpair> e;
    
    Eigen::eig(n, FmJac.data(), GmJac.data(), e);

    //if (family < 0 || family > e.size() - 1) return FIELD_ERROR;

    // TODO: Check also that lambda is not complex!

    double lambdam = e[family].r;
    current_diff_lambda = lambda_init_base_rarefaction - lambdam;
    
    return BISECTION_FUNCTION_OK;    
}

int CompositeCurve::sigma_minus_lambda_signal_event(const RealVector &where, double &sigma_minus_lambda, int *obj, int * /*not used*/){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    // The values associated with the rarefaction point.
    //
    double lambda_at_double_contact = composite_object->lambda_at_double_contact;

    RealVector rar_F_at_double_contact = composite_object->rar_F_at_double_contact; 
    RealVector rar_G_at_double_contact = composite_object->rar_G_at_double_contact;

    // Sigma. Extract the composite part of the field.
    //
    int two_n = where.size();
    int n = two_n/2;

    RealVector composite_point(n, n, where);

    JetMatrix cmp_F_jet(n), cmp_G_jet(n);
    composite_object->flux->jet(composite_point, cmp_F_jet, 0);
    composite_object->accum->jet(composite_point, cmp_G_jet, 0);

    double sigma = composite_object->shock->get_HugoniotContinuation()->sigma(rar_F_at_double_contact, 
                                                                              rar_G_at_double_contact, 
                                                                              cmp_F_jet.function(), 
                                                                              cmp_G_jet.function());

    sigma_minus_lambda = sigma - lambda_at_double_contact;

    return BISECTION_FUNCTION_OK; 
}

int CompositeCurve::sigma_minus_maxsigma_signal_event(const RealVector &where, double &sigma_minus_maxsigma, int *obj, int * /*not used*/){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    int two_n = where.size();
    int n = two_n/2;

    RealVector composite_point(n, n, where);

    JetMatrix cmp_F_jet(n), cmp_G_jet(n);
    composite_object->flux->jet(composite_point, cmp_F_jet, 0);
    composite_object->accum->jet(composite_point, cmp_G_jet, 0);

    // Sigma proper.
    //
    double sigma = composite_object->shock->get_HugoniotContinuation()->sigma(composite_object->rar_F_at_double_contact,
                                                                              composite_object->rar_G_at_double_contact,
                                                                              cmp_F_jet.function(), 
                                                                              cmp_G_jet.function());

    sigma_minus_maxsigma = sigma - composite_object->maxsigma;

    return BISECTION_FUNCTION_OK; 
}

int CompositeCurve::explicit_bifurcation_expression_signal_event(const RealVector &where, double &expression, int *obj, int * /*not used*/){
    CompositeCurve *composite_object = (CompositeCurve*)obj;

    int two_n = where.size();
    int n = two_n/2;

    RealVector composite_point(n, n, where);
    RealVector f = composite_object->explicit_bifurcation_curve->expressions(composite_point);

    expression = f(composite_object->index_of_explicit_bifurcation_expression);

    std::cout << "where           = " << where << std::endl;
    std::cout << "composite_point = " << composite_point << std::endl;
    std::cout << "f               = " << f << std::endl;
    std::cout << "index           = " << composite_object->index_of_explicit_bifurcation_expression << std::endl;
    std::cout << "expression      = " << expression << std::endl;

//    TestTools::pause();

    return BISECTION_FUNCTION_OK;
}

// TODO: As soon as possible have fp come from the outside, and be replaced with fq.
//
int CompositeCurve::transition_with_explicit_bifurcation(const ODE_Solver *odesolver, const RealVector &rarcmp_point, double init_time, RealVector &out, double &final_time){
//    class alpha_index {
//        public:
//            double alpha;
//            int index;

//            alpha_index() : alpha(0.0), index(0){}

//            alpha_index(double a, int i) : alpha(a), index(i){}

//            ~alpha_index(){}

//            bool operator<(const alpha_index &a){
//                return alpha < a.alpha;
//            }

//            alpha_index & operator=(const alpha_index &orig){
//                if (&orig != this){
//                    alpha = orig.alpha;
//                    index = orig.index;
//                }

//                return *this;
//            }
//    };

    int two_n = rarcmp_point.size();
    int n = two_n/2;

    RealVector fp = explicit_bifurcation_curve->expressions(RealVector(n, n, rarcmp_point));

    RealVector fq = explicit_bifurcation_curve->expressions(RealVector(n, n, out));

    std::vector<alpha_index> alphaindex;

    for (int k = 0; k < fq.size(); k++) {
        if (fp(k)*fq(k) < 0.0)  alphaindex.push_back(alpha_index(    -fq(k)/(fp(k) - fq(k)), k    ));
    }

    if (alphaindex.size() > 0){
        std::sort(alphaindex.rbegin(), alphaindex.rend());

//        {
//            std::stringstream ss;
//            for (int ii = 0; ii < alphaindex.size(); ii++) ss << "Alpha = " << alphaindex[ii].alpha << ", index = " << alphaindex[ii].index << std::endl;
//            TestTools::pause(ss);
//        }

        double bisection_epsilon = 1e-7;
        double c_t;
        RealVector p_c;

        index_of_explicit_bifurcation_expression = alphaindex[0].index;

        int info_bisection = Bisection::bisection_method(init_time,  rarcmp_point,
                                                         final_time, out,
                                                         bisection_epsilon, 
                                                         c_t, p_c,
                                                         &composite_field, (int*)this, (double*)0,
                                                         odesolver,
                                                         &explicit_bifurcation_expression_signal_event, (int*)this /*int *signal_event_object*/, (int*)0 /*int *signal_event_data*/);

        if (info_bisection == BISECTION_FUNCTION_OK){
            out = p_c;
            final_time = c_t;
        }

        return SECUNDARY_BIFURCATION_DETECTED;
    } 
    else {
        return SECUNDARY_BIFURCATION_NOT_DETECTED;
    }
}


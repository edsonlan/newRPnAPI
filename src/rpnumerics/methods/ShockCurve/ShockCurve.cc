#include "ShockCurve.h"

ShockCurve::ShockCurve(HugoniotContinuation *h){
    hc = h;

    boundary = hc->boundary();

    // TODO: This will change in the future: the physical domain and the computational domain
    //       need not be the same.
    //
    computational_domain = boundary;

    f = hc->flux();
    g = hc->accumulation();

    
}

ShockCurve::~ShockCurve(){
}

// This method is based on sigma_dot. If in a given interval there is an even number of changes of the sign of sigma_dot,
// this method will not detect them. Therefor, this method must be replaced by one using (sigma - lambda). The family will be needed.
//
void ShockCurve::find_system_for_sigma_equal_current_lambda(const RealVector &in, DoubleMatrix &nablaH, RealVector &H){
    int n = in.size();

    RealVector F(n), G(n);
    DoubleMatrix JF(n, n), JG(n, n);
//    double J2F[n][n][n], J2G[n][n][n];

//    f->fill_with_jet(n, in.components(), 2, F.components(), JF.data(), &J2F[0][0][0]);
//    g->fill_with_jet(n, in.components(), 2, G.components(), JG.data(), &J2G[0][0][0]);

    JetMatrix JM_F(n), JM_G(n);
    f->jet(WaveState(in), JM_F, 2);
    g->jet(WaveState(in), JM_G, 2);

    F = JM_F.function();
    G = JM_G.function();

    JF = JM_F.Jacobian();
    JG = JM_G.Jacobian();

    RealVector diff_F = F - ref.F;
    RealVector diff_G = G - ref.G;

    double inv_den = 1.0/(diff_G*diff_G);

    hc->jet_Hugoniot(F, JF, G, JG, H, nablaH);

    // Compute sigma
    double sigma = hc->sigma(F, G);

    // Find dsigma_du
    //
    RealVector dsigma_du = inv_den*(diff_G*JF + diff_F*JG - 2.0*sigma*diff_G*JG);
    //std::cout << "dsigma_du (vectorial) = " << dsigma_du << std::endl;

    //for (int j = 0; j < n; j++){
    //    dsigma_du(j) = 0.0;
//
//        // For all variables
//        for (int i = 0; i < n; i++){
//            dsigma_du(j) += JF(i, j)*diff_G(i) + JG(i, j)*diff_F(i) - 2.0*sigma*diff_G(i)*JG(i, j);
//        }
//
//        dsigma_du(j) *= inv_den;
//    }
    
//    std::cout << "dsigma_du (by components) = " << dsigma_du << std::endl;

    // Compute nablaH's last column and H's last element:
    DoubleMatrix A = JF - sigma*JG;

    H.resize(n);
    H(n - 1) = det(A);

//    nablaH.resize(n, n);
//    for (int i = 0; i < n; i++){
//        DoubleMatrix dA_dui(n, n);

//        for (int j = 0; j < n; j++){
//            for (int k = 0; k < n; k++){
////                dA_dui(j, k) = J2F[i][j][k] - dsigma_du(i)*JG(j, k) - sigma*J2G[i][j][k];
//                dA_dui(j, k) = J2F[j][i][k] - dsigma_du(i)*JG(j, k) - sigma*J2G[j][i][k];
//            }
//        }

//        nablaH(i, n - 1) = derivative_det(A, dA_dui);
//    }

    nablaH.resize(n, n);
    for (int i = 0; i < n; i++){
        DoubleMatrix dA_dui(n, n);

//        for (int j = 0; j < n; j++){
//            for (int k = 0; k < n; k++){
//                dA_dui(j, k) = J2F[i][j][k] - dsigma_du(i)*JG(j, k) - sigma*J2G[i][j][k];
//                dA_dui(j, k) = J2F[j][i][k] - dsigma_du(i)*JG(j, k) - sigma*J2G[j][i][k];
//            }
//        }

        //std::cout << "Size of A =\n" << dA_dui.rows() << std::endl;

        dA_dui = JM_F.extract_matrix_from_Hessian(i) - dsigma_du(i)*JG - sigma*JM_G.extract_matrix_from_Hessian(i);

        //std::cout << "After:  dA_dui =\n" << dA_dui << std::endl;

        nablaH(i, n - 1) = derivative_det(A, dA_dui);
    }

    return;
}

// This method is based on sigma_dot. If in a given interval there is an even number of changes of the sign of sigma_dot,
// this method will not detect them. Therefor, this method must be replaced by one using (sigma - lambda). The family will be needed.
//
int ShockCurve::find_point_for_sigma_equal_current_lambda(const RealVector &in, RealVector &out){
    RealVector iteration_point(in);  std::cout << "Iteration point = " << iteration_point << std::endl;

    int max_it = 20;
    int iterations = 0;

    bool found_point = false;

    DoubleMatrix nablaH;
    RealVector H, deviation;

    while (iterations < max_it && !found_point){
        find_system_for_sigma_equal_current_lambda(iteration_point, nablaH, H);

        std::cout << "Inside Newton. nablaH =\n" << nablaH << "H = " << H << std::endl;

        int info_solve = solve(transpose(nablaH), H, deviation);

        if (info_solve == REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR){
            std::cout << "Error in Newton. nablaH =\n" << nablaH << "H = " << H << std::endl;

            return REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR;
        }

        iteration_point = iteration_point - deviation; // was: - deviation

        // Verify that the point found by the Newton method is inside the
        // computational domain.
        //
        if (!computational_domain->inside(iteration_point)){
            std::cout << "Newton method for sigma equal current lambda: point " << iteration_point << " is off-limits!" << std::endl;
            return REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR;
        }

        if (norm(deviation) < 1e-8) found_point = true;

        iterations++;
    }

    if (found_point){
        std::cout << "*** Newton method for sigma equal current lambda converged!\n" << "Initial point = " << in << std::endl << std::endl << std::endl;

        out = iteration_point;
        return SHOCKCURVE_NEWTON_CONVERGED;
    }
    else {
        
        std::cout << "Newton method for sigma equal current lambda did not converge. Deviation = " << deviation << std::endl << std::endl << std::endl;

        return SHOCKCURVE_NEWTON_DID_NOT_CONVERGE;
    }
}

// This method computes nablaH and H for the Newton method in ShockCurve::find_point_for_sigma_equal_reference_lambda().
//

void ShockCurve::find_system(const RealVector &in, const RealVector &rarefaction_point, double lambda_ref, DoubleMatrix &nablaH, RealVector &H){

    ReferencePoint rar_with_stuff(rarefaction_point, f, g, 0);
    find_system(in, rar_with_stuff, lambda_ref, nablaH, H);

    return;
}


void ShockCurve::find_system(const RealVector &in, const ReferencePoint &reference_point, double lambda_ref, DoubleMatrix &nablaH, RealVector &H){
//    std::cout << "find system. in = " << in << std::endl;
/*
    int n = in.size();

    RealVector F(n), G(n);
    DoubleMatrix JF(n, n), JG(n, n);
    double J2F[n][n][n], J2G[n][n][n];

    f->fill_with_jet(n, in.components(), 2, F.components(), JF.data(), &J2F[0][0][0]);
    g->fill_with_jet(n, in.components(), 2, G.components(), JG.data(), &J2G[0][0][0]);
    
    
*/

    int n = in.size();

    RealVector F(n), G(n);
    DoubleMatrix JF(n, n), JG(n, n);

    JetMatrix JM_F(n), JM_G(n);
    f->jet(WaveState(in), JM_F, 2);
    g->jet(WaveState(in), JM_G, 2);

    F = JM_F.function();
    G = JM_G.function();

    JF = JM_F.Jacobian();
    JG = JM_G.Jacobian();

//    std::cout << "Inside find_system.\n" << "JM_F =\n" << JM_F << "JM_G = \n" << JM_G << std::endl;

//    std::cout << "Will invoke jet_Hugoniot" << std::endl;

    hc->jet_Hugoniot(F, JF, G, JG, H, nablaH);

    // Compute sigma
    double sigma = hc->sigma(F, G);

//    std::cout << "===>    Computed sigma: " << sigma << std::endl;
//    std::cout << "===>    Fref = " << ref.F << ", Gref = " << ref.G << std::endl;
//    std::cout << "===>      in = " << in << std::endl;
//    std::cout << "===>       F = " <<     F << ",    G = " << G << std::endl;
//    for (int i = 0; i < F.size(); i++) std::cout << "===>    Theoretical sigma = " << (F(i) - ref.F(i))/(G(i) - ref.G(i)) << ", component " << i << ". [F] = " << (F(i) - ref.F(i)) << ". [G] = " << (G(i) - ref.G(i)) << std::endl;

    // Print the derivative of sigma (delete later).
    {
        // By components.
        for (int i = 0; i < n; i++){
            double delta = .01;

            RealVector in_next(in);
            in_next(i) += delta;

            RealVector Fnext(n), Gnext(n);

            f->fill_with_jet(n, in_next.components(), 0, Fnext.components(), 0, 0);
            g->fill_with_jet(n, in_next.components(), 0, Gnext.components(), 0, 0);
            
            double sigma_next = (Fnext(i) - ref.F(i))/(Gnext(i) - ref.G(i));
            double sigma_now  = (F(i) - ref.F(i))/(G(i) - ref.G(i));

//            std::cout << "    sigma\'(" << i << ") = " << (sigma_next - sigma_now)/delta << std::endl;
        }
    }
    // Print the derivative of sigma (delete later).

    RealVector diff_F = F - ref.F;
    RealVector diff_G = G - ref.G;

    double inner_prod_FG = diff_F*diff_G;
    double inner_prod_GG = diff_G*diff_G;

    double inv_diff_G_norm_squared = 1.0/inner_prod_GG;
    double inv_diff_G_norm_squared_squared = inv_diff_G_norm_squared*inv_diff_G_norm_squared;

    // Prepare Newton proper
//    H(n - 1) = sigma - lambda_ref;
    H.resize(n);
    H(n - 1) = (sigma - lambda_ref)*(diff_G*diff_G);
//    H(n - 1) = (sigma - lambda_ref); // Perhaps the Newton method converges to a non-desired zero because we deal with a 3-order zero polynomial? This is supposed to fix it.

    // By rows, fill the last column.    
    nablaH.resize(n, n);
    for (int j = 0; j < n; j++){
        double num = 0.0;

        // For all variables
        for (int i = 0; i < n; i++){
            num += JF(i, j)*diff_G(i) + JG(i, j)*diff_F(i) - 2.0*lambda_ref*diff_G(i)*JG(i, j);
        }

        nablaH(j, n - 1) = num;
//        nablaH(j, n - 1) = num/(diff_G*diff_G); // Perhaps the Newton method converges to a non-desired zero because we deal with a 3-order zero polynomial? This is supposed to fix it.
                
    }

    return;
}

int ShockCurve::find_point_for_sigma_equal_reference_lambda(const RealVector &in, double lambda_ref, RealVector &out){
    RealVector iteration_point(in);
    hc->set_reference_point(ref);

//    std::cout << "Iteration point = " << iteration_point << ", reference point = " << ref.point << std::endl;

    int max_it = 20;
    int iterations = 0;

    bool found_point = false;

    DoubleMatrix nablaH;
    RealVector H, deviation;

    while (iterations < max_it && !found_point){
        find_system(iteration_point, ref, lambda_ref, nablaH, H);
        //find_system(iteration_point, lambda_ref, nablaH, H);

//        std::cout << "Inside Newton. nablaH =\n" << nablaH << "H = " << H << std::endl;

        int info_solve = solve(transpose(nablaH), H, deviation);

        if (info_solve == REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR){
            std::cout << "Error in Newton. nablaH =\n" << nablaH << "H = " << H << std::endl;

            return REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR;
        }

        iteration_point = iteration_point - deviation; // was: - deviation
//        std::cout << "Newton. Iteration = " << iterations << ", nablaH =\n" << nablaH << "det = " << det(nablaH) << ", deviation =\n   " << deviation << std::endl << "solution = " << iteration_point << std::endl;


        // Verify that the point found by the Newton method is inside the
        // computational domain.
        //
        if (!computational_domain->inside(iteration_point)){
            // Out must be filled even if iteration_point lies outside of the domain.
            out = iteration_point;

            std::cout << "Newton method within ShockCurve: point " << iteration_point << " is off-limits!" << std::endl;
            return SHOCKCURVE_NEWTON_OUTSIDE_DOMAIN;
        }

        if (norm(deviation) < 1e-8) found_point = true;

        iterations++;

           
    }

    std::cout << "*** Found_point = " << found_point << std::endl << std::endl << std::endl;

    if (found_point){
        out = iteration_point;
        return SHOCKCURVE_NEWTON_CONVERGED;
    }
    else {
        
        std::cout << "Newton did not converge. Iterations: " << iterations << std::endl;

        return SHOCKCURVE_NEWTON_DID_NOT_CONVERGE;
    }
}

// TODO: After checking the Bethe-Wendroff Theorem, lambda_minus_sigma can be replaced by sigma_dot.
//       Limitations: If sigma_dot vanishes an even number of times, this method will not detect them.
//       Therefor, using lambda - sigma could be advantageous.
//
void ShockCurve::find_alphas_for_characteristic_shocks(const RealVector &previous_lambda_minus_sigma, 
                                                       const RealVector &lambda_minus_sigma, 
                                                       int what_family_to_use, int continue_after_transition, 
                                                       std::vector<double> &alpha,
                                                       std::vector<int> &corresponding_family){
    alpha.clear();
    corresponding_family.clear();

    // Find all alphas...
    //
    if (what_family_to_use == USE_ALL_FAMILIES){
        for (int i = 0; i < previous_lambda_minus_sigma.size(); i++){
            if (previous_lambda_minus_sigma(i)*lambda_minus_sigma(i) < 0.0){
                alpha.push_back(lambda_minus_sigma(i)/(previous_lambda_minus_sigma(i) - lambda_minus_sigma(i)));
                corresponding_family.push_back(i);
            }
        }

    // Alternatively, only for the case where a zero of sigma_dot is to be found,
    // because the BW Theorem guarantees that sigma_dot = 0 for one and only one family.

//        //
//        int pos = 0;
//        bool found = false;

//        while (pos < previous_lambda_minus_sigma.size() && !found){
//            if (previous_lambda_minus_sigma(pos)*lambda_minus_sigma(pos) < 0.0){
//                alpha.push_back(lambda_minus_sigma(pos)/(previous_lambda_minus_sigma(pos) - lambda_minus_sigma(pos)));
//                found = true;
//            }

//            pos++;
//        }
//        //

        // Sort in reverse order.
        //
        std::sort(alpha.rbegin(), alpha.rend());

        // If the curve is to stop after the transition, only the point closest to the previous point will be
        // used.
        //
        if (continue_after_transition == STOP_AFTER_TRANSITION){
            if (alpha.size() > 0) alpha.resize(1);
        }
    }
    // ...or only one alpha
    else if (what_family_to_use == USE_CURRENT_FAMILY){
        if (previous_lambda_minus_sigma(current_family)*lambda_minus_sigma(current_family) < 0.0){
            alpha.push_back(lambda_minus_sigma(current_family)/(previous_lambda_minus_sigma(current_family) - lambda_minus_sigma(current_family)));
            corresponding_family.push_back(current_family);
        }
    }

    return;
}

// FIRST TEST.
// Check if sigma_dot changes sign.
//
// This method is allowed to modify Hugoniot_intersection and Hugoniot_direction.
//
// By the Bethe-Wendroff Theorem, sigma_dot = 0 iff there exists a family such that 
//
//     sigma = lambda(family).
//
// This may affect (slightly) the jet_Hugoniot methods.
//
int ShockCurve::local_speed_equality(const RealVector &previous_lambda_minus_sigma, const RealVector &previous_point, 
                                     const RealVector &lambda_minus_sigma,          const RealVector &candidate_point,
                                     int what_family_to_use,
                                     std::vector<TransitionPointStructure> &transition_points){
//    std::cout << "local_speed_equality: previous = " << previous_point << ", candidate = " << candidate_point << std::endl;

    std::vector<double> alpha;
    std::vector<int>    corresponding_family;

    find_alphas_for_characteristic_shocks(previous_lambda_minus_sigma, 
                                          lambda_minus_sigma, 
                                          what_family_to_use, continue_after_transition, 
                                          alpha, corresponding_family);

    for (int i = 0; i < alpha.size(); i++){
        // Dan hated this one here. alpha's meaning should be redefined to (1.0 - alpha), so it is bound to the last point and not to the previous point.
        RealVector approximate_transition = (1.0 - alpha[i])*candidate_point + alpha[i]*previous_point;

        RealVector accurate_transition;

        int info = find_point_for_sigma_equal_current_lambda(approximate_transition, accurate_transition);

        if (info == SHOCKCURVE_NEWTON_CONVERGED){
            transition_points.push_back(TransitionPointStructure(alpha[i], accurate_transition, corresponding_family[i], true));
            std::cout << "local_speed_equality." << std::endl;
            for (int i = 0; i < transition_points.size(); i++) std::cout << "    local = " << transition_points[i].local << std::endl;
        }
        //else SHOCKCURVE_NEWTON_DID_NOT_CONVERGE;
    }

//    //Hugoniot_intersection = (1.0 - alpha)*Hugoniot_intersection + alpha*shockcurve[shockcurve.size() - 1];
//    candidate_point = (1.0 - alpha)*candidate_point + alpha*previous_point; /*** HERE ***/

//    int info_find_point_for_sigma_equal_current_lambda = find_point_for_sigma_equal_current_lambda(RealVector(candidate_point), candidate_point);

//    // Update Hugoniot_direction
//    hc->find_continuation_direction(candidate_point, RealVector(Hugoniot_direction), Hugoniot_direction);
//    transition_current.push_back(candidate_point);

    return CONTINUE_AFTER_TRANSITION;
}

// SECOND TEST.
//
// Check if sigma == lambda_ref for any family.
//
// This method is allowed to modify Hugoniot_intersection and Hugoniot_direction.
//
//int ShockCurve::reference_speed(RealVector &Hugoniot_intersection, RealVector &Hugoniot_direction){
//    double alpha_sigma_minus_lambda = 1.0;
//    bool found_sigma_minus_lambda = false;
//    double lambda_ref_family_to_be_used;

//    for (int i = 0; i < sigma_minus_lambda.size(); i++){
//        JetMatrix sigmaj(1);
//        hc->sigma_jet(Hugoniot_intersection, (const RealVector *)0, 0, sigmaj);
//        double sigma_Hugoniot = sigmaj.get(0);
//    
//        double lambda_ref_family = ref.e[i].r;

//        if ((sigma_Hugoniot - lambda_ref_family)*sigma_minus_lambda[i] < 0.0){
//            double alpha_temp = (sigma_Hugoniot - lambda_ref_family)/((sigma_Hugoniot - lambda_ref_family) - sigma_minus_lambda[i]);

//            if (alpha_temp < alpha_sigma_minus_lambda) alpha_sigma_minus_lambda = alpha_temp;
//            lambda_ref_family_to_be_used = lambda_ref_family;

//            found_sigma_minus_lambda = true;

//        }

//        sigma_minus_lambda[i] = sigma_Hugoniot - lambda_ref_family; // Lambda_ref
//    }

//    if (found_sigma_minus_lambda){
//        Hugoniot_intersection = alpha_sigma_minus_lambda*shockcurve[shockcurve.size() - 1] + (1.0 - alpha_sigma_minus_lambda)*Hugoniot_intersection;

//        // Use Newton
//        RealVector temp_h(Hugoniot_intersection);

//        int info_find_point_for_sigma_equal_reference_lambda = find_point_for_sigma_equal_reference_lambda(temp_h, lambda_ref_family_to_be_used, Hugoniot_intersection);

//        // Update Hugoniot_direction
//        hc->find_continuation_direction(Hugoniot_intersection, RealVector(Hugoniot_direction), Hugoniot_direction);

//        transition_reference.push_back(Hugoniot_intersection);

//        return STOP_AFTER_TRANSITION;
//    }
//    else CONTINUE_AFTER_TRANSITION;
//}

int ShockCurve::reference_speed_equality(const RealVector &previous_lambdaref_minus_sigma, const RealVector &previous_point, 
                                         const RealVector &lambdaref_minus_sigma,          const RealVector &candidate_point,
                                         int what_family_to_use,
                                         std::vector<TransitionPointStructure> &transition_points){
    std::vector<double> alpha;
    std::vector<int>    corresponding_family;

    find_alphas_for_characteristic_shocks(previous_lambdaref_minus_sigma, 
                                          lambdaref_minus_sigma, 
                                          what_family_to_use, continue_after_transition, 
                                          alpha, corresponding_family);

    for (int i = 0; i < alpha.size(); i++){
        // Dan hated this one here. alpha's meaning should be redefined to (1.0 - alpha), so it is bound to the last point and not to the previous point.
        RealVector approximate_transition = (1.0 - alpha[i])*candidate_point + alpha[i]*previous_point;

        RealVector accurate_transition;

        int info = find_point_for_sigma_equal_reference_lambda(approximate_transition, ref.e[corresponding_family[i]].r, accurate_transition);

        if (info == SHOCKCURVE_NEWTON_CONVERGED){
            transition_points.push_back(TransitionPointStructure(alpha[i], accurate_transition, corresponding_family[i], false));
            std::cout << "reference_speed_equality." << std::endl;
            for (int i = 0; i < transition_points.size(); i++) std::cout << "    local = " << transition_points[i].local << std::endl;
        }
        //else SHOCKCURVE_NEWTON_DID_NOT_CONVERGE;
    }

    return CONTINUE_AFTER_TRANSITION;
}

void ShockCurve::add_point(Curve &c, const RealVector &p){
    c.curve.push_back(p);

    // The corresponding point of All points 
    // in the shock curve is the initial point of the initial rarefaction curve.
    c.back_pointer.push_back(0); 

    int n = p.size();
    JetMatrix Fjet(n), Gjet(n);
    
    f->jet(p, Fjet, 1);
    g->jet(p, Gjet, 1);

//    std::cout << "Fjet = \n" << Fjet << std::endl << std::endl;
//    std::cout << "Gjet = \n" << Gjet << std::endl << std::endl;

//    std::cout << "Fjet.Jacobian() = \n" << Fjet.Jacobian() << std::endl << std::endl;
//    std::cout << "Gjet.Jacobian() = \n" << Gjet.Jacobian() << std::endl << std::endl;

    if (Fjet.Jacobian().cols() == 0){
        std::cout << "At shock. Fjet.Jacobian().cols() = 0. ABORTING." << std::endl;
        std::cout << "    p = " << p << std::endl;
        exit(0);
    }

    // Eigenvalues.
    //
    std::vector<eigenpair> e;

    Eigen::eig(n, Fjet.Jacobian().data(), Gjet.Jacobian().data(), e);

    RealVector eigenvalues(e.size());
    for (int i = 0; i < e.size(); i++) eigenvalues(i) = e[i].r;

    c.eigenvalues.push_back(eigenvalues);

    // Speed. If within a certain distance of the reference point, use the eigenvalue.
    //
    if (c.curve.size() == 1) c.speed.push_back(e[current_family].r);
    else                     c.speed.push_back(hc->sigma(Fjet.function(), Gjet.function()));

//    std::cout << "Adding shock point. p = " << p << ", speed = " << c.speed.back() << std::endl;

    return;
}

// TODO: Generalize the interruption to say if (lambda - sigma) increases or decreases, at each family used.
// TODO: Write classifier based on the interruption functions. Possible answers: "Lax shock of some family", "Classification in Alexei-Marchesin".
// TODO: Write classifier without using the interruption functions for non-connected branches.
//
//
int ShockCurve::call_interruption_functions(const RealVector &previous_lambdaref_minus_sigma, const RealVector &previous_lambda_minus_sigma, const RealVector &previous_point, 
                                            const RealVector &lambdaref_minus_sigma,          const RealVector &lambda_minus_sigma,          const RealVector &candidate_point,
                                            int what_family_to_use,
                                            int after_transition,
                                            int left_subtype, int right_subtype,
                                            Curve &curve,          // To be updated by this method, adding the transition points.
                                            std::vector<int> &transition_current_index,
                                            std::vector<int> &transition_current_family,
                                            int &transition_current_found,
                                            std::vector<int> &transition_reference_index,
                                            std::vector<int> &transition_reference_family,
                                            int &transition_reference_found){

//    // TODO: These clears are WRONG!
//    transition_current_index.clear();
//    transition_current_family.clear();
//    transition_reference_index.clear();
//    transition_reference_family.clear();

    std::vector<TransitionPointStructure> transition_points;

    if (right_subtype != DONT_CHECK_EQUALITY_AT_RIGHT){
        local_speed_equality(previous_lambda_minus_sigma, previous_point,
                             lambda_minus_sigma,          candidate_point,
                             what_family_to_use,
                             transition_points);

    }

    if (left_subtype != DONT_CHECK_EQUALITY_AT_LEFT){
        reference_speed_equality(previous_lambdaref_minus_sigma, previous_point,
                                 lambdaref_minus_sigma,          candidate_point,
                                 what_family_to_use,
                                 transition_points);
    }

    // All transition points are sorted here.
    //
    std::sort(transition_points.begin(), transition_points.end());

    // Redistribute the transition points according to their being local or with respect to the reference point.
    //
    transition_current_found = transition_reference_found = TRANSITION_NOT_FOUND;

    for (int i = 0; i < transition_points.size(); i++){
        int n = curve.curve.size();

        if (transition_points[i].local == true){
            transition_current_index.push_back(n);
            transition_current_family.push_back(transition_points[i].family);
            transition_current_found = TRANSITION_CURRENT_FOUND;

            std::cout << "Added one to current." << std::endl;
        }
        else{
            transition_reference_index.push_back(n);
            transition_reference_family.push_back(transition_points[i].family);
            transition_reference_found = TRANSITION_REFERENCE_FOUND;

            std::cout << "Added one to reference." << std::endl;
        }

        add_point(curve, transition_points[i].point);
        //curve.curve.push_back(transition_points[i].point);
    }

    return (transition_points.size() > 0) ? TRANSITION_FOUND : TRANSITION_NOT_FOUND;
}

int ShockCurve::curve_engine(const ReferencePoint &r, const RealVector &in, const RealVector &initial_direction, int family, 
                             int type, int left_subtype, int right_subtype,
                             int what_family_to_use,
                             int after_transition,
                             void *linobj, double (*linear_function)(void *o, const RealVector &p),
                             Curve &shockcurve, 
                             std::vector<int> &transition_current_index,
                             std::vector<int> &transition_current_family,
                             std::vector<int> &transition_reference_index,
                             std::vector<int> &transition_reference_family,
                             int &shock_stopped_because,
                             int &edge){

    shockcurve.clear();
    shockcurve.type = SHOCK_CURVE;
    shockcurve.family = family;
    shockcurve.back_curve_index = 0; // This can be changed by the WaveCurveFactory. In principle it's ok.

    transition_current_index.clear();
    transition_current_family.clear();
    transition_reference_index.clear();
    transition_reference_family.clear();

    // For find_alphas_for_characteristic_shocks().
    current_family = family;

    ref = r;
    hc->set_reference_point(r);

    RealVector previous_point = in;

    add_point(shockcurve, in); // See relation with lambda.

    RealVector previous_direction = initial_direction;
    normalize(previous_direction);

    // Initialize the intersection of the rarefaction and a line.
    //
    double old_linear_function_value = 0.0;
    if (linear_function != 0) old_linear_function_value = (*linear_function)(linobj, previous_point);
    double linear_function_value = old_linear_function_value;

    double step_size = hc->default_step_size();
    int number_of_steps_with_unchanged_size = 0;
    int step_size_increased = 0;

    const Boundary *b = hc->boundary();

//    double lambda = r.e[family].r;
//    double sigma_minus_lambda = 0.0;

    double sigma_dot_old = 0.0; // Use d_lambda!
    double sigma_old;

    int n = r.e.size();

    std::vector<double> sigma_minus_lambda(r.e.size());
    for (int i = 0; i < sigma_minus_lambda.size(); i++) sigma_minus_lambda[i] = 0.0;

    // Lambda - sigma. Size depends on what_family_to_use.
    //
    int size_lambda_minus_sigma;

    if (what_family_to_use == USE_CURRENT_FAMILY) size_lambda_minus_sigma = 1;
    else                                          size_lambda_minus_sigma = n;

    RealVector previous_lambda_minus_sigma(size_lambda_minus_sigma);
    RealVector lambda_minus_sigma(size_lambda_minus_sigma);

    // Lambda_ref - sigma.
    //
    RealVector previous_lambdaref_minus_sigma(n);
    RealVector lambdaref_minus_sigma(n), candidate_point;

    // To be used by HugoniotContinuation::curve_point() to determine if the step_size can be changed or not.
    //
    double previous_sigma_between_points = 0.0;

    while (true){
        std::cout << previous_point << std::endl;

        RealVector Hugoniot_intersection;
        double sigma_between_points;
        RealVector Hugoniot_direction;

//        std::cout << "Shock. Prev. point = " << previous_point << ", prev. dir. = " << previous_direction << std::endl;

        int info_curve_point = hc->curve_point(previous_point, previous_sigma_between_points,
                                               previous_direction, 
                                               step_size_increased, step_size, number_of_steps_with_unchanged_size, 
                                               Hugoniot_intersection, sigma_between_points,
                                               Hugoniot_direction);

//        std::cout << "       Curr. point = " << Hugoniot_intersection << ", curr. dir. = " << Hugoniot_direction << std::endl;
//        std::cout << "       sigma_between_points = " << sigma_between_points << std::endl;

        // Update sigma_between_points
        previous_sigma_between_points = sigma_between_points;

        // Verify if the new point lies within the domain or not.
        // TODO: Why is it that this is not working correctly for the Liquid?
        //

//        previous_point     = Hugoniot_intersection;
//        previous_direction = Hugoniot_direction;

        double sigma, sigma_dot;
        hc->jet_sigma(Hugoniot_intersection, Hugoniot_direction, sigma, sigma_dot);

        candidate_point = Hugoniot_intersection;

//        // Stop criteria below. If at least one of them is met, the computation of the curve is stopped.
//        bool add_last_point_computed = true;
//        bool find_next_point         = true;

        // Update previous_lambdaref_minus_sigma & lambdaref_minus_sigma.
        //
        previous_lambdaref_minus_sigma = lambdaref_minus_sigma;
            
        for (int i = 0; i < n; i++) lambdaref_minus_sigma(i) = r.e[i].r - sigma;

        // Update previous_lambda_minus_sigma & lambda_minus_sigma.
        //
        previous_lambda_minus_sigma = lambda_minus_sigma;
         
        std::vector<double> lambda;
        Eigen::fill_eigenvalues(f, g, candidate_point, lambda);

        if (what_family_to_use == USE_CURRENT_FAMILY) lambda_minus_sigma(0) = lambda[family] - sigma;
        else {
            lambda_minus_sigma.resize(lambda.size()); // TODO: See if this can be moved upwards. Is it possible that the number of eigenvalues changes for the same model?
            for (int i = 0; i < lambda.size(); i++) lambda_minus_sigma(i) = lambda[i] - sigma;
        }
 
        if (shockcurve.curve.size() > 1){
            // Specific stop criteria. The caller MUST set use_interruption_functions.
            // TODO: I think it is best to pass use_interruption_functions in the signature.
            //
            if (use_interruption_functions != DONT_USE_INTERRUPTION_FUNCTIONS) {
                int transition_current_found, transition_reference_found;

                int info_interruption = call_interruption_functions(previous_lambdaref_minus_sigma, previous_lambda_minus_sigma, previous_point,
                                                                    lambdaref_minus_sigma,          lambda_minus_sigma,          candidate_point,
                                                                    what_family_to_use,
                                                                    after_transition,
                                                                    left_subtype, right_subtype,
                                                                    shockcurve,
                                                                    transition_current_index,
                                                                    transition_current_family,
                                                                    transition_current_found,
                                                                    transition_reference_index,
                                                                    transition_reference_family,
                                                                    transition_reference_found);
                //std::cout << "After call_interruptions. Current transitions: " << transition_current_index.size() << ", reference transitions: " << transition_reference_index.size() << std::endl;

                // TODO: Depending on the type of transition the curve must stop if after_transition == STOP_AFTER_TRANSITION.
                //       Now it is simply stopping if any transition is found.
                //
                if (info_interruption == TRANSITION_FOUND && after_transition == STOP_AFTER_TRANSITION){
                    if (transition_current_found == TRANSITION_CURRENT_FOUND){
                        shock_stopped_because = SHOCK_RIGHT_CHARACTERISTIC_AT_FAMILY;
                    }

                    else if (transition_reference_found == TRANSITION_REFERENCE_FOUND){
                        shock_stopped_because = SHOCK_LEFT_CHARACTERISTIC_AT_FAMILY;
                    }

                    shockcurve.final_direction = shockcurve.curve.back() - shockcurve.curve[shockcurve.curve.size() - 2];
                    normalize(shockcurve.final_direction);

                    shockcurve.last_point = shockcurve.curve.back();

                    return SHOCKCURVE_OK; //TODO: What to return here?
                }
            }

        } //  if (shockcurve.size() > 1)

        if (shockcurve.curve.size() > 0){
            RealVector r;
            int info_intersect = b->intersection(previous_point, Hugoniot_intersection, r, edge);
//            std::cout << "Intersection: r = " << r << std::endl;

            if (linear_function != 0){
                linear_function_value = (*linear_function)(linobj, Hugoniot_intersection);

                if (linear_function_value*old_linear_function_value < 0.0){
                    double alpha = linear_function_value/(linear_function_value - old_linear_function_value);

                    RealVector point_on_line = alpha*previous_point + (1.0 - alpha)*Hugoniot_intersection;

                    add_point(shockcurve, point_on_line);

                    shock_stopped_because = SHOCK_REACHED_LINE;
                    shockcurve.reason_to_stop = SHOCK_REACHED_LINE;

                    shockcurve.last_point = point_on_line;
                    shockcurve.final_direction = Hugoniot_intersection - previous_point; //rarcurve.curve.back() - rarcurve.curve[rarcurve.curve.size() - 2];
                    normalize(shockcurve.final_direction);

                    return SHOCKCURVE_OK;
                }
                else old_linear_function_value = linear_function_value;
            }

            // Both points are inside: carry on.
            if (info_intersect == BOUNDARY_INTERSECTION_BOTH_INSIDE){
//                add_last_point_computed = true;
                //find_next_point         = true;

                add_point(shockcurve, Hugoniot_intersection);
            }
            // Both outside (this really should not happen).
//            else if (info_intersect == BOUNDARY_INTERSECTION_BOTH_OUTSIDE){
//                return HUGONIOTCONTINUATION_CURVE_ERROR;
//            }
            // New point outside: the curve reached the domain's boundary.
            else {  // info_intersect == BOUNDARY_INTERSECTION_FOUND
//                Hugoniot_intersection = r;

//                add_last_point_computed = true;
//                find_next_point         = false;                

                std::cout << "Shock. Curve leaves domain:" << std::endl;
                std::cout << "    previous_point = " << previous_point << std::endl;
                std::cout << "             point = " << Hugoniot_intersection << std::endl;
                std::cout << "              info = " << info_intersect << std::endl;

                add_point(shockcurve, r);

                shock_stopped_because = SHOCK_REACHED_BOUNDARY;
                return SHOCKCURVE_OK;

                //return HUGONIOTCONTINUATION_CURVE_OK;
            }
        }

        // NEW LINE BELOW, update Hugoniot_direction using Hugoniot_intersection.
        //
        hc->find_continuation_direction(Hugoniot_intersection, RealVector(Hugoniot_direction), Hugoniot_direction);

        sigma_old     = sigma;
        sigma_dot_old = sigma_dot;

        previous_point     = Hugoniot_intersection;
        previous_direction = Hugoniot_direction;
    }
}

int ShockCurve::curve_engine_from_boundary(RarefactionCurve *rarefactioncurve, int side, int increase, const ReferencePoint &r, const RealVector &in, int family, 
                                 int type, int left_subtype, int right_subtype,
                                 int what_family_to_use,
                                 int after_transition,
                                 Curve &shockcurve, 
                                 std::vector<int> &transition_current_index,
                                 std::vector<int> &transition_current_family,
                                 std::vector<int> &transition_reference_index,
                                 std::vector<int> &transition_reference_family, 
                                 int &shock_stopped_because,
                                 int &edge){

    // Points to the interior of the domain from side s.
    //
    RealVector to_interior = boundary->side_transverse_interior(in, side);

    // Initialize.
    //
    RealVector initial_direction;
    double dd;

    int info_initialize = rarefactioncurve->initialize(in, family, increase, initial_direction, dd);

    if (info_initialize == RAREFACTION_INIT_ERROR) return SHOCKCURVE_ERROR;

    // Follow the direction opposite to the rarefaction's.
    //
    if (initial_direction*to_interior < 0.0){
        int info = curve_engine(r, in, initial_direction, family, 
                                type, left_subtype, right_subtype,
                                what_family_to_use,
                                after_transition,
                                shockcurve, 
                                transition_current_index,
                                transition_current_family,
                                transition_reference_index,
                                transition_reference_family, 
                                shock_stopped_because,
                                edge);

        return info;
    }
    else {
        return SHOCKCURVE_ERROR;
    }
}


int ShockCurve::curve(const ReferencePoint &ref, 
                      int type, int left_subtype, int right_subtype,
                      int what_family_to_use, int after_transition,
                      std::vector<ShockCurvePoints> &shockcurve_with_transitions){                
//                      std::vector< std::vector<RealVector> > &curve, 
//                      std::vector< std::vector<RealVector> > &transition_current,
//                      std::vector< std::vector<RealVector> > &transition_reference,
//                      std::vector<ShockCurvePoints> shockcurve_with_transitions){

//    curve.clear();
    shockcurve_with_transitions.clear();

    int n = ref.point.size();

    std::cout << "Here" << std::endl;

//    // If the initial point lies near the coincidence curve, abort.
//    // TODO: These lines below need to be improved.
//    //
//    if (reference_point_near_coincidence()){ // TODO: This value must be fine-tuned.

////        initial_step_size = ...;
//        return HUGONIOTCONTINUATION_NEAR_COINCIDENCE_CURVE;
//    }

    for (int family = 0; family < ref.e.size(); family++){
        // Find the eigenvector of this family.
        RealVector r(n);
        for (int i = 0; i < n; i++) r(i) = ref.e[family].vrr[i];

        RealVector final_direction;
        int edge;

        Curve shockcurve;
        std::vector<int> transition_current_index;
        std::vector<int> transition_current_family;
        std::vector<int> transition_reference_index;
        std::vector<int> transition_reference_family;
        int shock_stopped_because;

//        curve_engine(ref, ref.point, r, family, type, subtype, temp, temp_transition_current, temp_transition_reference, edge);

        curve_engine(ref, ref.point, r, family, type, left_subtype, right_subtype, 
                     what_family_to_use,
                     after_transition,
                     shockcurve, 
                     transition_current_index,
                     transition_current_family,
                     transition_reference_index,
                     transition_reference_family,  
                     shock_stopped_because,
                     edge);

        shockcurve_with_transitions.push_back(ShockCurvePoints(shockcurve.curve, 
                                                               transition_current_index,   transition_current_family,
                                                               transition_reference_index, transition_reference_family,
                                                               family, ref));

//        if (temp.size() > 0) curve.push_back(temp);
//        if (temp_transition_current.size() > 0)   transition_current.push_back(temp_transition_current);
//        if (temp_transition_reference.size() > 0) transition_reference.push_back(temp_transition_reference);

//        temp.clear();
//        temp_transition_current.clear();
//        temp_transition_reference.clear();

//        curve_engine(ref, ref.point, -r, family, type, subtype, temp, temp_transition_current, temp_transition_reference, edge);

//        if (temp.size() > 0) curve.push_back(temp);
//        if (temp_transition_current.size() > 0)   transition_current.push_back(temp_transition_current);
//        if (temp_transition_reference.size() > 0) transition_reference.push_back(temp_transition_reference);

        shockcurve.clear();
        transition_current_index.clear();
        transition_current_family.clear();
        transition_reference_index.clear();
        transition_reference_family.clear();  

        curve_engine(ref, ref.point, -r, family, type, left_subtype, right_subtype,
                     what_family_to_use,
                     after_transition,
                     shockcurve, 
                     transition_current_index,
                     transition_current_family,
                     transition_reference_index,
                     transition_reference_family,  
                     shock_stopped_because,
                     edge);

        shockcurve_with_transitions.push_back(ShockCurvePoints(shockcurve.curve, 
                                                               transition_current_index,   transition_current_family,
                                                               transition_reference_index, transition_reference_family,
                                                               family, ref));

    }

    return HUGONIOTCONTINUATION_CURVE_OK;
}

// TODO: THIS METHOD IS INCOMPLETE!!!
//
// There are two complementary ways to code this method:
//
//     a) The curve that contains the Bethe-Wendroff points was created using this ShockCurve's flux, accumulation, etc., while the curve_containing_the_extension was created with different flux, accumulation, etc.,
//     b) The curve_containing_the_extension was created using this ShockCurve's flux, accumulation, etc., while the curve that contains the Bethe-Wendroff points was created with different flux, accumulation, etc.
//
// The first alternative is in use here. The method could also be coded as a static method, but that would break the style.
//
//                                               ==================================
//                                               NOTE ON HOW THE RESULTS ARE STORED
//                                               ==================================
//
// The extension of a Bethe-Wendroff point may be formed by several points. Since this method extends several Bethe-Wendroff points,
// it would be more natural to store the results in
//
//     std::vector<std::vector<RealVector> > &extension,
//
// that is, a vector of points per Bethe-Wendroff point. On the downside, if no extension points are found for a particular Bethe-Wendroff
// point, the corresponding std::vector<RealVector> would be empty. To avoid this situation the results are stored as a pair of vectors:
// of points and indices:
//
//     std::vector<RealVector> &extension      and     std::vector<int> &corresponding_Bethe_Wendroff.
//
// In "extension" all the extended points will be stored, while in "corresponding_Bethe_Wendroff" will be stored the indices of the corresponding
// Bethe-Wendroff point that was extended.
//
//
//                                               ==============================
//                                               NOTE ON THE RESULTS THEMSELVES
//                                               ==============================
//
// There seems to be a problem with the results (the extension), which is found where it should not be. Perhaps it is a problem with the updating
// of the sigmas. THIS SHOULD BE TAKEN CARE OF! To complete this method we are waiting for a correct/precise definition of "extension of a Bethe-Wendroff point".
//
void ShockCurve::Bethe_Wendroff_extension(const Curve &curve, const std::vector<int> &transition_current_index, const std::vector<int> &transition_current_family, const Curve &curve_where_extension_is_to_be_found, 
                                          std::vector<RealVector> &curve_with_extension, std::vector<int> &corresponding_Bethe_Wendroff){

    // Create a copy.
    //
    curve_with_extension.clear();
    for (int i = 0; i < curve_where_extension_is_to_be_found.curve.size(); i++){
        curve_with_extension.push_back(curve_where_extension_is_to_be_found.curve[i]);
    }

    corresponding_Bethe_Wendroff.clear();

    const FluxFunction *ext_f;
    const AccumulationFunction *ext_g;

    // This may change in the future, so far all the computations will happen on the same domain.
    ext_f = f;
    ext_g = g;
    HugoniotContinuation *ext_hc = hc;

    // TODO: This ShockCurve may come from the outside. I don't like this solution.
    //
    ShockCurve ext_sc(ext_hc);

    for (int i = 0; i < transition_current_index.size(); i++){
        int n = curve.curve[transition_current_index[i]].size();

        JetMatrix F_J(n), G_J(n);
        f->jet(curve.curve[transition_current_index[i]], F_J, 0);
        g->jet(curve.curve[transition_current_index[i]], G_J, 0);

        RealVector F = F_J.function();
        RealVector G = G_J.function();

        double Bethe_Wendroff_speed = curve.speed[transition_current_index[i]];

        JetMatrix ext_F_J(n), ext_G_J(n);
        ext_f->jet(curve_where_extension_is_to_be_found.curve[0], ext_F_J, 0);
        ext_g->jet(curve_where_extension_is_to_be_found.curve[0], ext_G_J, 0);

        double prev_sigma = ext_hc->sigma(F, G, ext_F_J.function(), ext_G_J.function());
//        std::cout << curve.curve[transition_current_index[i]] << ", Bethe-Wendroff sigma = " << Bethe_Wendroff_speed << std::endl;

        for (int j = 1; j < curve_where_extension_is_to_be_found.curve.size(); j++){
            // If curve and curve_to_be_extended are the same, avoid extending 
            if ((&curve != &curve_where_extension_is_to_be_found) || (transition_current_index[i] != j)){
                ext_f->jet(curve_where_extension_is_to_be_found.curve[j], ext_F_J, 0);
                ext_g->jet(curve_where_extension_is_to_be_found.curve[j], ext_G_J, 0);

                double sigma = ext_hc->sigma(F, G, ext_F_J.function(), ext_G_J.function());

                if ((sigma - Bethe_Wendroff_speed)*(prev_sigma - Bethe_Wendroff_speed) < 0.0){
                    // Compute the point. TODO: 
                    //
                    RealVector extension_point = .5*(curve_where_extension_is_to_be_found.curve[j] + curve_where_extension_is_to_be_found.curve[j - 1]);

                    curve_with_extension.insert(curve_with_extension.begin() + j, extension_point);
                    corresponding_Bethe_Wendroff.push_back(transition_current_index[i]);

//                    std::cout << "    **** Extension found:" << std::endl;
                }

//                std::cout << "    point = " << curve_to_be_extended.curve[j] << std::endl;
//                std::cout << "    prev_sigma = " << prev_sigma << ", sigma = " << sigma << std::endl << std::endl;

                prev_sigma = sigma;
            }
        }
    }

    return;
}

// TODO: Test this method.
void ShockCurve::Bethe_Wendroff_extension(const Curve &curve, 
                                          const std::vector<int> &transition_current_index, 
                                          std::vector<ExtensionPoint> &curve_plus_extensions){
    curve_plus_extensions.clear();

    // Make a copy of the curve that is being extended.
    //
    for (int i = 0; i < curve.curve.size(); i++) curve_plus_extensions.push_back(ExtensionPoint(curve.curve[i], -1));

    for (int i = 0; i < transition_current_index.size(); i++){
        double Bethe_Wendroff_speed = curve.speed[transition_current_index[i]];

        double prev_sigma = curve.speed[0];

        int tci = transition_current_index[i];

        for (int j = 1; j < curve.curve.size(); j++){
            double sigma = curve.speed[j];

            // Avoid extending a point onto itself.
            //
            if (j != tci){
                if ((sigma - Bethe_Wendroff_speed)*(prev_sigma - Bethe_Wendroff_speed) < 0.0){
                    // Compute the point. TODO: For starters, use a linear interpolation.
                    //
                    double alpha = (Bethe_Wendroff_speed - sigma)/(prev_sigma - sigma);
                    RealVector approximate_extension = alpha*curve.curve[j - 1] + (1.0 - alpha)*curve.curve[j];
                    RealVector accurate_extension(approximate_extension);

//                    // TODO: For some reason the method below is failing to compute the accurate extension.
//                    int info_extension = find_point_for_sigma_equal_reference_lambda(approximate_extension, Bethe_Wendroff_speed, accurate_extension);
//                    if (info_extension != SHOCKCURVE_NEWTON_CONVERGED) accurate_extension = approximate_extension;

                    ExtensionPoint extension_point(accurate_extension, tci);

                    curve_plus_extensions.insert(curve_plus_extensions.begin() + j, extension_point);
                }
            }

            prev_sigma = sigma;
        }
    }

    // TODO: The output must change: the method must return a Curve.

    return;
}


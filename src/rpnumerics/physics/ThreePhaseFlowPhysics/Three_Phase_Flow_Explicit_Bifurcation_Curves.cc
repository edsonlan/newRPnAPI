#include "Three_Phase_Flow_Explicit_Bifurcation_Curves.h"

/* Ok */
Three_Phase_Flow_Explicit_Bifurcation_Curves::Three_Phase_Flow_Explicit_Bifurcation_Curves(ThreePhaseFlowFluxFunction *f) : Explicit_Bifurcation_Curves(), flux(f) {
}

/* Ok */
Three_Phase_Flow_Explicit_Bifurcation_Curves::~Three_Phase_Flow_Explicit_Bifurcation_Curves(){
}

/* Ok */
RealVector Three_Phase_Flow_Explicit_Bifurcation_Curves::sec_bif_correspondence(int side_opposite_vertex, const RealVector &point, const RealVector &mu){
    int n = point.size();

    double D0 = 0.0;
    for (int i = 0; i < n; i++) D0 += point(i)*point(i)/mu(i);

    // Temporal
    double sg = 1.0 - point(0) - point(1);
    D0 += sg*sg/mu(2);

    RealVector correspondence(3);

    if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SW_ZERO){
        double temp = 1.0 + mu(2)/mu(1);

        correspondence(1) = (point(1)/mu(0))/(D0 + (temp*temp/mu(0) + temp/mu(1))*(point(1)*point(1)));
        correspondence(2) = correspondence(1)*mu(2)/mu(1);
        correspondence(0) = 1.0 - correspondence(1) - correspondence(2);
    }
    else if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SO_ZERO){
        double temp = 1.0 + mu(2)/mu(0);

        correspondence(0) = (point(0)/mu(1))/(D0 + (temp*temp/mu(1) + temp/mu(0))*(point(0)*point(0)));
        correspondence(2) = correspondence(0)*mu(2)/mu(0);
        correspondence(1) = 1.0 - correspondence(0) - correspondence(2);
    }
    else {
        double temp = 1.0 + mu(1)/mu(0);

        correspondence(0) = (point(0)/mu(2))/(D0 + (temp*temp/mu(2) + temp/mu(0))*(point(0)*point(0)));
        correspondence(1) = correspondence(0)*mu(1)/mu(0);
        correspondence(2) = 1.0 - correspondence(0) - correspondence(1);
    } 

    return correspondence;
}

/* Ok */
void Three_Phase_Flow_Explicit_Bifurcation_Curves::vertex_and_side(int side_opposite_vertex, const RealVector &mu, RealVector &vertex, RealVector &point_on_side){
    vertex.resize(2);
    point_on_side.resize(2);

    double muw = mu(0);
    double muo = mu(1);
    double mug = mu(2);

    if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SW_ZERO){
        vertex(0) = 1.0;
        vertex(1) = 0.0;

        point_on_side(0) = 0.0;
        point_on_side(1) = muo/(mug + muo); // Equation of line from            sg/mug = so/muo
                                            //                       (1 - sw - so)/mug = so/muo, or, (1 - sw - so)*muo - so*mug = 0.
    }
    else if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SO_ZERO){
        vertex(0) = 0.0;
        vertex(1) = 1.0;

        point_on_side(0) = muw/(muw + mug); // Equation of line from            sw/muw = sg/mug,
                                            //                                  sw/muw = (1 - sw - so)/mug, or, sw*mug - (1 - sw - so)*muw = 0.
        point_on_side(1) = 0.0;
    }
    else {
        vertex(0) = 0.0;
        vertex(1) = 0.0;

        point_on_side(0) = muw/(muw + muo); // Equation of line from            sw/muw = so/muo, or, sw*muo - so*muw = 0.
        point_on_side(1) = muo/(muw + muo);
    }  

    return;
}

/* Ok */
// Concrete method declared in the parent class.
//
int Three_Phase_Flow_Explicit_Bifurcation_Curves::sec_bif_correspondence(int side_opposite_vertex, int nos, 
                                                                         std::vector<RealVector> &point, 
                                                                         std::vector<RealVector> &correspondent_point,
                                                                         int &model_specific_error_code){

    int info_permeability_parameters = check_permeability_parameters(model_specific_error_code);
    if (info_permeability_parameters == PERMEABILITY_PARAMETERS_ERROR){
        return EXPLICIT_BIFURCATION_CODE_ERROR;
    }

    int info_flux_parameters = check_flux_parameters(model_specific_error_code);
    if (info_flux_parameters == FLUX_PARAMETERS_ERROR){
        return EXPLICIT_BIFURCATION_CODE_ERROR;
    }

    // Extract viscosities.
    //

    RealVector mu(3);
    mu(0) = flux->muw()->value(); // muw
    mu(2) = flux->mug()->value(); // mug
    mu(1) = flux->muo()->value(); // muo

    // Endpoints
    RealVector vertex(2), point_on_side(2);

    vertex_and_side(side_opposite_vertex, mu, vertex, point_on_side); 

    // Just in case...
    //
    nos = std::max(2, nos);

    RealVector delta = (point_on_side - vertex)/((double)nos - 1);

    for (int i = 0; i < nos; i++){
        RealVector p = vertex + delta*((double)i);

        point.push_back(p);
        correspondent_point.push_back(sec_bif_correspondence(side_opposite_vertex, p, mu));
    }

    // Temporary kludge
    Utilities::points_to_segments(point);
    Utilities::points_to_segments(correspondent_point);
    
    return EXPLICIT_BIFURCATION_CODE_OK;
}

/* Ok */
RealVector Three_Phase_Flow_Explicit_Bifurcation_Curves::expressions(const RealVector &point){
    // muw, muo, mug

    double muw = flux->muw()->value();
    double mug = flux->mug()->value();
    double muo = flux->muo()->value();

    double sw = point(0);
    double so = point(1);

    RealVector eq(3);

    eq(0) = so*mug - (1.0 - sw - so)*muo;
    eq(1) = sw*mug - (1.0 - sw - so)*muw;
    eq(2) = so*muw - sw*muo;

    return eq;
}

/* Ok */
int Three_Phase_Flow_Explicit_Bifurcation_Curves::region(const RealVector &expressions){
    double eq_SW = expressions(0);
    double eq_SO = expressions(1);
    double eq_SG = expressions(2);

    if      (eq_SW < 0.0 && eq_SO < 0.0 && eq_SG < 0.0) return REGION_WM_OM_GM;
    else if (eq_SW < 0.0 && eq_SO > 0.0 && eq_SG < 0.0) return REGION_WM_OP_GM;
    else if (eq_SW > 0.0 && eq_SO > 0.0 && eq_SG < 0.0) return REGION_WP_OP_GM;
    else if (eq_SW > 0.0 && eq_SO > 0.0 && eq_SG > 0.0) return REGION_WP_OP_GP;
    else if (eq_SW > 0.0 && eq_SO < 0.0 && eq_SG > 0.0) return REGION_WP_OM_GP;
    else if (eq_SW < 0.0 && eq_SO < 0.0 && eq_SG > 0.0) return REGION_WM_OM_GP;
}

//int Three_Phase_Flow_Explicit_Bifurcation_Curves::find_point_for_sigma_equal_current_lambda(const RealVector &in, RealVector &out){
//    RealVector iteration_point(in);

//    int max_it = 20;
//    int iterations = 0;

//    bool found_point = false;

//    DoubleMatrix nablaH;
//    RealVector H, deviation;

//    while (iterations < max_it && !found_point){
//        shock->find_system_for_sigma_equal_current_lambda(iteration_point, nablaH, H);

//        std::cout << "Inside Newton. nablaH =\n" << nablaH << "H = " << H << std::endl;

//        int info_solve = solve(transpose(nablaH), H, deviation);

//        if (info_solve == REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR){
//            std::cout << "Error in Newton. nablaH =\n" << nablaH << "H = " << H << std::endl;

//            return REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR;
//        }

//        iteration_point = iteration_point - deviation; // was: - deviation

//        // Verify that the point found by the Newton method is inside the
//        // computational domain.
//        //
//        if (!computational_domain->inside(iteration_point)){
//            std::cout << "Newton method for sigma equal current lambda: point " << iteration_point << " is off-limits!" << std::endl;
//            return REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR;
//        }

//        if (norm(deviation) < 1e-8) found_point = true;

//        iterations++;
//    }

//    if (found_point){
//        std::cout << "*** Newton method for sigma equal current lambda converged!\n" << "Initial point = " << in << std::endl << std::endl << std::endl;

//        out = iteration_point;
//        return SHOCKCURVE_NEWTON_CONVERGED;
//    }
//    else {
//        
//        std::cout << "Newton method for sigma equal current lambda did not converge. Deviation = " << deviation << std::endl << std::endl << std::endl;

//        return SHOCKCURVE_NEWTON_DID_NOT_CONVERGE;
//    }
//}

/* Ok */
void Three_Phase_Flow_Explicit_Bifurcation_Curves::subdivide_curve_in_regions(const Curve &curve, 
                                                                              Curve &classified_curve,
                                                                              std::vector<int> &Lax_transition_index,
                                                                              std::vector<int> &region_transition_index){
    classified_curve.clear();
    region_transition_index.clear();

    if (curve.curve.size() < 2) return;

    RealVector fp, fq;

    // Initialization.
    //
    fp = expressions(curve.curve[0]);
    int n = fp.size();
    classified_curve.curve.push_back(curve.curve[0]);

    // Iterate. 
    //
    for (int i = 1; i < curve.curve.size(); i++){
        fq = expressions(curve.curve[i]);

        std::vector<double> alpha;
        
        for (int j = 0; j < n; j++) if (fp(j)*fq(j) < 0.0) alpha.push_back(-fq(j)/(fp(j) - fq(j)));

        if (alpha.size() > 0){
            std::sort(alpha.begin(), alpha.end());

            // Store the transition points.
            //
            for (int j = 0; j < alpha.size(); j++){
                classified_curve.curve.push_back(   alpha[j]*curve.curve[i - 1] + (1.0 - alpha[j])*curve.curve[i]   );

//                region_transition_index.push_back(i - 1 + j);
//                region_transition_index.push_back(i + j);
                region_transition_index.push_back(classified_curve.curve.size() - 1);
            }

            // Update Lax_transition_index.
            //
            for (int j = 0; j < Lax_transition_index.size(); j++){
                if (i > Lax_transition_index[j]) Lax_transition_index[j] += alpha.size();
            }
        } 

        // Store the endpoint of this segment.
        //
        classified_curve.curve.push_back(curve.curve[i]);

        // Update.
        //
        fp = fq;
    }

    return;
}

void Three_Phase_Flow_Explicit_Bifurcation_Curves::subdivide_segmented_curve_in_regions(const std::vector<RealVector>&curve, 
                                                                                        std::vector<RealVector> &classified_curve,
                                                                                        std::vector<int> &Lax_transition_index,
                                                                                        std::vector<int> &region_transition_index){
    classified_curve.clear();
    region_transition_index.clear();

    if (curve.size() < 2) return;

    for (int i = 0; i < curve.size()/2; i++){
        Curve temp_classified_curve;
        std::vector<int> temp_Lax_transition_index;
        std::vector<int> temp_region_transition_index;

        Curve c;
        for (int j = 0; j < 2; j++) c.curve.push_back(curve[2*i + j]);

        subdivide_curve_in_regions(c, temp_classified_curve, temp_Lax_transition_index, temp_region_transition_index);

        // This line is NOT interchangeable with the one below. region_transition_index MUST be updated
        // before classified_curve.size() changes.
        //
        // Other possibility (that I think could allow swapping the lines): 
        //
        //     for (int j = 0; j < temp_region_transition_index.size(); j++) region_transition_index.push_back(2*i + region_transition_index.size() + temp_region_transition_index[j]);
        //
        for (int j = 0; j < temp_region_transition_index.size(); j++) region_transition_index.push_back(classified_curve.size() + temp_region_transition_index[j]);

        // Convert to segments.
        //
        Utilities::points_to_segments(temp_classified_curve.curve);

        for (int j = 0; j < temp_classified_curve.curve.size(); j++){
            classified_curve.push_back(temp_classified_curve.curve[j]);
        }
    }

    return;
}

int Three_Phase_Flow_Explicit_Bifurcation_Curves::check_permeability_parameters(int &model_specific_error_code){  
    return PERMEABILITY_PARAMETERS_OK;
}


int Three_Phase_Flow_Explicit_Bifurcation_Curves::check_flux_parameters(int &model_specific_error_code){
    return FLUX_PARAMETERS_OK;
}


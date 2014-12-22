#include "Stone_Explicit_Bifurcation_Curves.h"


/* OK */
Stone_Explicit_Bifurcation_Curves::Stone_Explicit_Bifurcation_Curves(StoneFluxFunction *ff) : Three_Phase_Flow_Explicit_Bifurcation_Curves((ThreePhaseFlowFluxFunction*)ff), f(ff){
}

/* OK */
Stone_Explicit_Bifurcation_Curves::~Stone_Explicit_Bifurcation_Curves(){
}

/* OK */
int Stone_Explicit_Bifurcation_Curves::check_permeability_parameters(int &model_specific_error_code){
    double expw = f->permeability()->expw()->value();
    double expg = f->permeability()->expg()->value();
    double expo = f->permeability()->expo()->value();

    double expow = f->permeability()->expow()->value();
    double expog = f->permeability()->expog()->value();

    double cnw = f->permeability()->cnw()->value();
    double cng = f->permeability()->cng()->value();
    double cno = f->permeability()->cno()->value();

    double lw = f->permeability()->lw()->value();
    double lg = f->permeability()->lg()->value();

    double low = f->permeability()->low()->value();
    double log = f->permeability()->log()->value();

    double epsl = f->permeability()->epsl()->value();

    // Only compute the curves in this case.
    //
    if (expw  == 2.0 && expg  == 2.0 && expo == 2.0 &&
        expow == 2.0 && expog == 2.0 &&
        cnw   == 0.0 && cng   == 0.0 && cno == 0.0 &&
        low   == 0.0 && log   == 0.0 &&
        epsl  == 0.0 ) return PERMEABILITY_PARAMETERS_OK;
    else return PERMEABILITY_PARAMETERS_ERROR;
}

/* OK */
int Stone_Explicit_Bifurcation_Curves::check_flux_parameters(int &model_specific_error_code){
    double grw = f->grw()->value();
    double grg = f->grg()->value();
    double gro = f->gro()->value();

    // Only compute the curves in this case.
    //
    if (!(grw == grg && grw == gro)){
        return FLUX_PARAMETERS_ERROR;
    }
    else return FLUX_PARAMETERS_OK;
}

//void Stone_Explicit_Bifurcation_Curves::vertex_and_side(int side_opposite_vertex, const RealVector &mu, RealVector &vertex, RealVector &point_on_side){
//    vertex.resize(2);
//    point_on_side.resize(2);

//    double muw = mu(0);
//    double muo = mu(1);
//    double mug = mu(2);

//    if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SW_ZERO){
//        vertex(0) = 1.0;
//        vertex(1) = 0.0;

//        point_on_side(0) = 0.0;
//        point_on_side(1) = muo/(mug + muo); // Equation of line from            sg/mug = so/muo
//                                            //                       (1 - sw - so)/mug = so/muo, or, (1 - sw - so)*muo - so*mug = 0.
//    }
//    else if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SO_ZERO){
//        vertex(0) = 0.0;
//        vertex(1) = 1.0;

//        point_on_side(0) = muw/(muw + mug); // Equation of line from            sw/muw = sg/mug,
//                                            //                                  sw/muw = (1 - sw - so)/mug, or, sw*mug - (1 - sw - so)*muw = 0.
//        point_on_side(1) = 0.0;
//    }
//    else {
//        vertex(0) = 0.0;
//        vertex(1) = 0.0;

//        point_on_side(0) = muw/(muw + muo); // Equation of line from            sw/muw = so/muo, or, sw*muo - so*muw = 0.
//        point_on_side(1) = muo/(muw + muo);
//    }  

//    return;
//}

/* OK */
void Stone_Explicit_Bifurcation_Curves::line(const RealVector &p, const RealVector &q, int nos, std::vector<RealVector> &v){
    v.clear();

    RealVector delta = (q - p)/((double)nos);

    // Add first point (p).
    //
    v.push_back(p);

    // Not starting a 0!!!
    //
    for (int i = 1; i <= nos; i++){
        RealVector temp = p + delta*((double)i);

        // Add two copies of the new point...
        //
        v.push_back(temp);
        v.push_back(temp);
    }

    // ...and remove the second copy of the last point.
    //
    v.pop_back();

    return;
}

/* OK */
void Stone_Explicit_Bifurcation_Curves::expl_sec_bif_crv(int side_opposite_vertex, int nos, 
                                                         std::vector<RealVector> &vertex_to_umbilic, 
                                                         std::vector<RealVector> &umbilic_to_side){
    vertex_to_umbilic.clear();
    umbilic_to_side.clear();

    double expw = f->permeability()->expw()->value();
    double expg = f->permeability()->expg()->value();
    double expo = f->permeability()->expo()->value();

    double expow = f->permeability()->expow()->value();
    double expog = f->permeability()->expog()->value();

    double cnw = f->permeability()->cnw()->value();
    double cng = f->permeability()->cng()->value();
    double cno = f->permeability()->cno()->value();

    double lw = f->permeability()->lw()->value();
    double lg = f->permeability()->lg()->value();

    double low = f->permeability()->low()->value();
    double log = f->permeability()->log()->value();

    double epsl = f->permeability()->epsl()->value();

    // Only compute the curves in this case.
    //
    if (expw  == 2.0 && expg  == 2.0 && expo == 2.0 &&
        expow == 2.0 && expog == 2.0 &&
        cnw   == 0.0 && cng   == 0.0 && cno == 0.0 &&
        low   == 0.0 && log   == 0.0 &&
        epsl  == 0.0 ){

        // muw, muo, mug
        double muw = f->muw()->value();
        double mug = f->mug()->value();
        double muo = f->muo()->value();

        RealVector mu(3);
        mu(0) = muw;
        mu(1) = muo;
        mu(2) = mug;

        double sum_mu = muw + muo + mug;

        // Umbilic point.
        //
        RealVector umbp(2);
        umbp(0) = muw/sum_mu;
        umbp(1) = muo/sum_mu;

        // Just in case...
        //
        nos = std::max(2, nos);

        // Endpoints
        RealVector vertex(2), point_on_side(2);

        vertex_and_side(side_opposite_vertex, mu, vertex, point_on_side);

        line(vertex, umbp, nos, vertex_to_umbilic);
        line(umbp, point_on_side, nos, umbilic_to_side);
    }

    return;
}

//RealVector Stone_Explicit_Bifurcation_Curves::expression(const RealVector &point){
//    // muw, muo, mug
//    RealVector flux_params = f->fluxParams().params();

//    double muw = flux_params(3);
//    double mug = flux_params(4);
//    double muo = flux_params(5);

//    double sw = point(0);
//    double so = point(1);

//    RealVector eq(3);

//    eq(0) = so*mug - (1.0 - sw - so)*muo;
//    eq(1) = sw*mug - (1.0 - sw - so)*muw;
//    eq(2) = so*muw - sw*muo;

//    return eq;
//}

//int Stone_Explicit_Bifurcation_Curves::region(const RealVector &equations){
//    double eq_SW = equations(0);
//    double eq_SO = equations(1);
//    double eq_SG = equations(2);

//    if      (eq_SW < 0.0 && eq_SO < 0.0 && eq_SG < 0.0) return REGION_WM_OM_GM;
//    else if (eq_SW < 0.0 && eq_SO > 0.0 && eq_SG < 0.0) return REGION_WM_OP_GM;
//    else if (eq_SW > 0.0 && eq_SO > 0.0 && eq_SG < 0.0) return REGION_WP_OP_GM;
//    else if (eq_SW > 0.0 && eq_SO > 0.0 && eq_SG > 0.0) return REGION_WP_OP_GP;
//    else if (eq_SW > 0.0 && eq_SO < 0.0 && eq_SG > 0.0) return REGION_WP_OM_GP;
//    else if (eq_SW < 0.0 && eq_SO < 0.0 && eq_SG > 0.0) return REGION_WM_OM_GP;
//}

//int Stone_Explicit_Bifurcation_Curves::cross_sec_bif(const RealVector &previous_point, const RealVector &point, RealVector &crossing_point, int &region){
//    // muw, muo, mug
//    RealVector flux_params = f->fluxParams().params();

//    double muw = flux_params(3);
//    double mug = flux_params(4);
//    double muo = flux_params(5);

//    double sw_prev = previous_point(0);
//    double so_prev = previous_point(1);

//    double sw      = point(0);
//    double so      = point(1);

//    double eq_W = so*mug - (1.0 - sw - so)*muo;
//    double eq_O = sw*mug - (1.0 - sw - so)*muw;
//    double eq_G = so*muw - sw*muo;

//    if      (eq_W < 0.0 && eq_O < 0.0 && eq_G < 0.0) region = REGION_WM_OM_GM; // ---
//    else if (eq_W < 0.0 && eq_O > 0.0 && eq_G < 0.0) region = REGION_WM_OP_GM; // -+-
//    else if (eq_W > 0.0 && eq_O > 0.0 && eq_G < 0.0) region = REGION_WP_OP_GM; // ++-
//    else if (eq_W > 0.0 && eq_O > 0.0 && eq_G > 0.0) region = REGION_WP_OP_GP; // +++
//    else if (eq_W > 0.0 && eq_O < 0.0 && eq_G > 0.0) region = REGION_WP_OM_GP; // +-+
//    else if (eq_W < 0.0 && eq_O < 0.0 && eq_G > 0.0) region = REGION_WM_OM_GP; // --+
//    
//    return 1;
//}

//RealVector Stone_Explicit_Bifurcation_Curves::sec_bif_correspondence(int side_opposite_vertex, const RealVector &point, const RealVector &mu){
//    int n = point.size();

//    double D0 = 0.0;
//    for (int i = 0; i < n; i++) D0 += point(i)*point(i)/mu(i);

//    // Temporal
//    double sg = 1.0 - point(0) - point(1);
//    D0 += sg*sg/mu(2);

//    RealVector correspondence(3);

//    if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SW_ZERO){
//        double temp = 1.0 + mu(2)/mu(1);

////        correspondence(1) = (point(1)/mu(1))/(D0 + (temp*temp/mu(0) + temp/mu(1))*(point(1)*point(1)));
//        correspondence(1) = (point(1)/mu(0))/(D0 + (temp*temp/mu(0) + temp/mu(1))*(point(1)*point(1)));
//        correspondence(2) = correspondence(1)*mu(2)/mu(1);
//        correspondence(0) = 1.0 - correspondence(1) - correspondence(2);
//    }
//    else if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SO_ZERO){
//        double temp = 1.0 + mu(2)/mu(0);

//        correspondence(0) = (point(0)/mu(1))/(D0 + (temp*temp/mu(1) + temp/mu(0))*(point(0)*point(0)));
//        correspondence(2) = correspondence(0)*mu(2)/mu(0);
//        correspondence(1) = 1.0 - correspondence(0) - correspondence(2);
//    }
//    else {
//        double temp = 1.0 + mu(1)/mu(0);

//        correspondence(0) = (point(0)/mu(2))/(D0 + (temp*temp/mu(2) + temp/mu(0))*(point(0)*point(0)));
//        correspondence(1) = correspondence(0)*mu(1)/mu(0);
//        correspondence(2) = 1.0 - correspondence(0) - correspondence(1);
//    } 

//    return correspondence;
//}

//int Stone_Explicit_Bifurcation_Curves::sec_bif_correspondence(int side_opposite_vertex, int nos, 
//                                                              std::vector<RealVector> &point, 
//                                                              std::vector<RealVector> &correspondent_point,
//                                                              int &model_specific_error_code){
//    point.clear();
//    correspondent_point.clear();

//    RealVector permeability_params = f->perm().params().params();

//    double expw = permeability_params(0);
//    double expg = permeability_params(1);
//    double expo = permeability_params(2);

//    double expow = permeability_params(3);
//    double expog = permeability_params(4);

//    double cnw = permeability_params(5);
//    double cng = permeability_params(6);
//    double cno = permeability_params(7);

//    double lw = permeability_params(8);
//    double lg = permeability_params(9);

//    double low = permeability_params(10);
//    double log = permeability_params(11);

//    double epsl = permeability_params(12);

//    // Only compute the curves in this case.
//    //
//    if (expw  == 2.0 && expg  == 2.0 && expo == 2.0 &&
//        expow == 2.0 && expog == 2.0 &&
//        cnw   == 0.0 && cng   == 0.0 && cno == 0.0 &&
//        low   == 0.0 && log   == 0.0 &&
//        epsl  == 0.0 ){

//        // muw, muo, mug
//        RealVector flux_params = f->fluxParams().params();

//        double grw = flux_params(0);
//        double grg = flux_params(1);
//        double gro = flux_params(2);

//        if (!(grw == grg && grw == gro)){
//            int &model_specific_error_code = 
//            return;
//        }

//        RealVector mu(3);
//        mu(0) = flux_params(3);
//        mu(2) = flux_params(4);
//        mu(1) = flux_params(5);

//        std::cout << "mu = " << mu << std::endl;

//        // Endpoints
//        RealVector vertex(2), point_on_side(2);

//        vertex_and_side(side_opposite_vertex, mu, vertex, point_on_side); 

//        // Just in case...
//        //
//        nos = std::max(2, nos);

//        RealVector delta = (point_on_side - vertex)/((double)nos - 1);

//        for (int i = 0; i < nos; i++){
//            RealVector p = vertex + delta*((double)i);

//            point.push_back(p);
//            correspondent_point.push_back(sec_bif_correspondence(side_opposite_vertex, p, mu));
//        }
//    }
//    else {
//        model_specific_error_code = INADEQUATE_PERMEABILITY_PARAMETERS;
//        return EXPLICIT_BIFURCATION_CODE_ERROR;
//    }

//    return EXPLICIT_BIFURCATION_CODE_OK;
//}

//int Stone_Explicit_Bifurcation_Curves::check_flux_parameters(int &model_specific_error_code){
//    RealVector flux_params = f->fluxParams().params();

//    double grw = flux_params(0);
//    double grg = flux_params(1);
//    double gro = flux_params(2);

//    // Only compute the curves in this case.
//    //
//    if (!(grw == grg && grw == gro)){
//        return FLUX_PARAMETERS_ERROR;
//    }
//    else return FLUX_PARAMETERS_OK;
//}


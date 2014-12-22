#include "WaveCurveFactory.h"

// Find the intesection between two segments: p1-p2 and q1-q2, store the answer in r.
// If there is no intersection, return false (and r is useless), otherwise return true.
//
bool WaveCurveFactory::segment_intersection(double *p1, double *p2, double *q1, double *q2, double *r){
    double alpha, beta;

    double A[2][2], b[2];
    for (int i = 0; i < 2; i++){
        A[i][0] = p1[i] - p2[i];
        A[i][1] = q2[i] - q1[i];

        b[i]    = q2[i] - p2[i];
    }

    double delta = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (fabs(delta) < 1e-10) {
        return false;
    }

    alpha = (b[0]*A[1][1] - b[1]*A[0][1])/delta;
    beta  = (b[1]*A[0][0] - b[0]*A[1][0])/delta;

    for (int i = 0; i < 2; i++) r[i] = .5*(alpha*p1[i] + (1.0 - alpha)*p2[i] + beta*q1[i] + (1.0 - beta)*q2[i]);

    return (alpha >= 0.0 && alpha <= 1.0) && (beta >= 0.0 && beta <= 1.0);
}

WaveCurveFactory::WaveCurveFactory(const AccumulationFunction *gg, const FluxFunction *ff, const Boundary *bb, const ODE_Solver *o,
                                   RarefactionCurve *r, ShockCurve *s, CompositeCurve *c){
    g = gg;
    f = ff;
    b = bb;

    odesolver = o;

    rarefactioncurve = r;
    shockcurve       = s;
    compositecurve   = c;

    type.resize(3);
    type[0] = std::string("Rarefaction Curve");
    type[1] = std::string("Composite Curve");
    type[2] = std::string("Shock Curve");
}

WaveCurveFactory::~WaveCurveFactory(){
}

int WaveCurveFactory::Liu_half_wavecurve(const ReferencePoint &ref, 
                                         const RealVector &initial_point, 
                                         int initial_family, 
                                         int increase, 
                                         int initial_curve, 
                                         const RealVector &initial_direction, 
                                         void *linobj, double (*linear_function)(void *o, const RealVector &p),
                                         WaveCurve &hwc,
                                         int &wavecurve_stopped_because, 
                                         int &edge){

    // Kludge to solve the fact that the shockspeed at the reference point is being returned as NaN, 
    // which messes up with the Riemann Profile.
    //
    bool is_first = true;

    // This changes if the rarefaction reaches a coincidence.
    //
    int family = initial_family;

    int future_curve = initial_curve;
    RealVector future_curve_initial_point(initial_point);
    RealVector future_curve_initial_direction(initial_direction);

    int future_composite_type;

    std::vector<int> rarefaction_list;
    std::vector<int> last_point_in_rarefaction;

    while (true){
        if (hwc.wavecurve.size() > 3) return WAVECURVE_OK;

        int current_curve = future_curve;
        RealVector current_curve_initial_point(future_curve_initial_point);
        RealVector current_curve_initial_direction(future_curve_initial_direction);

        std::cout << "WaveCurveFactory, top of while cycle.\n    current_curve_initial_point = " <<  current_curve_initial_point << "\n    current_curve_initial_direction = " << current_curve_initial_direction << std::endl;

        if (current_curve == RAREFACTION_CURVE){
            std::cout << "WaveCurveFactory: entering Rarefaction." << std::endl;

            double deltaxi = 1e-3;//3e-4; // Was: 1e-3
            std::vector<RealVector> inflection_point;
            Curve rarcurve;

            int rar_stopped_because;
            RealVector final_direction;

            int info_rar = rarefactioncurve->curve(current_curve_initial_point,
                                                   family,
                                                   increase,
                                                   RAREFACTION,
                                                   RAREFACTION_INITIALIZE, //RAREFACTION_DONT_INITIALIZE,
                                                   &current_curve_initial_direction,
                                                   odesolver,
                                                   deltaxi,
                                                   linobj, linear_function,
                                                   rarcurve,
                                                   inflection_point,
                                                   final_direction,
                                                   rar_stopped_because,
                                                   edge);

            // Update the back pointers.
            //
            rarcurve.back_curve_index = hwc.wavecurve.size() - 1;

            is_first = false;

            if (rarcurve.curve.size() > 0){
                rarcurve.back_curve_pointer.resize(rarcurve.curve.size());
                rarcurve.back_curve_pointer[0] = hwc.wavecurve.size() - 1;
                for (int i = 1; i < rarcurve.curve.size(); i++) rarcurve.back_curve_pointer[i] = hwc.wavecurve.size();
                
                if (hwc.wavecurve.size() > 0) rarcurve.back_pointer[0] = hwc.wavecurve.back().curve.size() - 1;
            }

            // Store regardless of what happened during the computation.
            //
            hwc.wavecurve.push_back(rarcurve);

//            // Remove this!             
//            return WAVECURVE_OK;

            if (info_rar == RAREFACTION_OK){
                // The rarefaction is its own reference curve.
                //
                rarcurve.back_curve_index = hwc.wavecurve.size();
                rarcurve.reference_point = ref;

                future_curve_initial_point = rarcurve.last_point;
                future_curve_initial_direction = rarcurve.final_direction;

                if (rar_stopped_because == RAREFACTION_REACHED_INFLECTION){
                    future_curve = COMPOSITE_CURVE;
                    future_composite_type = COMPOSITE_BEGINS_AT_INFLECTION;

                    // Add this rarefaction to the stack of rarefactions that future composites will try to exhaust.
                    //
                    if (rarcurve.curve.size() > 0){
                        rarefaction_list.push_back(hwc.wavecurve.size() - 1);
                        last_point_in_rarefaction.push_back(rarcurve.curve.size() - 1);
                    }
                    else {
                        std::cout << "WaveCurveFactory: the rarefaction claims it finished ok, but returned zero points!" << std::endl;
                    }

                    std::cout << "RarefactionCurve stopped at inflection. List of rarefactions:" << std::endl;
                    for (int i = 0; i < rarefaction_list.size(); i++) std::cout << "    Rar = " << rarefaction_list[i] << ", index at rar = " << last_point_in_rarefaction[i] << std::endl;
                }
                else if (rar_stopped_because == RAREFACTION_COMPLEX_EIGENVALUE_AT_FAMILY){
                    wavecurve_stopped_because = WAVECURVE_COMPLEX_EIGENVALUE_AT_FAMILY;

                    std::cout << "WaveCurveFactory: exiting at Rarefaction (complex eigenvalue)." << std::endl;

                    return WAVECURVE_OK;
                }
                else if (rar_stopped_because == RAREFACTION_REACHED_COINCIDENCE_CURVE){
                    // Helmut-like only. In principle, it should continue as a rarefaction with the next family.
                    //
                    wavecurve_stopped_because = WAVECURVE_REACHED_COINCIDENCE_CURVE;

                    std::cout << "WaveCurveFactory: exiting at Rarefaction (coincidence curve)." << std::endl;

                    return WAVECURVE_OK;
                }
                else if (rar_stopped_because == RAREFACTION_REACHED_BOUNDARY){
                    wavecurve_stopped_because = WAVECURVE_REACHED_BOUNDARY;

                    // In multidomains this will change.

                    // There could be non-connected Composites associated with this rarefaction.
                    // In that case, something should be here.

                    std::cout << "WaveCurveFactory: exiting at Rarefaction (reached boundary)." << std::endl;

                    return WAVECURVE_OK;
                }
                else if (rar_stopped_because == RAREFACTION_REACHED_COINCIDENCE_CURVE){
                    future_curve = RAREFACTION_CURVE;
                    if (increase == RAREFACTION_SPEED_SHOULD_INCREASE) family = family + 1;
                    else                                               family = family - 1;

                    if (family >= initial_point.size() - 1 || family < 0){
                        return WAVECURVE_ERROR;
                    }

                    // Add this rarefaction to the stack of rarefactions that future composites will try to exhaust.
                    //
                    if (rarcurve.curve.size() > 0){
                        rarefaction_list.push_back(hwc.wavecurve.size() - 1);
                        last_point_in_rarefaction.push_back(rarcurve.curve.size() - 1);
                    }
                    else {
                        std::cout << "WaveCurveFactory: the rarefaction claims it finished ok, but returned zero points!" << std::endl;
                    }
                }
                else if (rar_stopped_because == RAREFACTION_REACHED_LINE){
                    wavecurve_stopped_because = WAVECURVE_REACHED_LINE;

                    std::cout << "WaveCurveFactory reached the line." << std::endl;

                    return WAVECURVE_OK;
                }
                // TODO: add Panters' case. lambda = 0!!!
//                else{
//                }
                else {
                    // This point should not be reached: it indicates that the rarefaction ended in
                    // a way that the WaveCurveFactory cannot handle. Warn the user and abort.

                    return WAVECURVE_ERROR;
                }
            }
            else { // if (info_rar == RAREFACTION_ERROR)
                return WAVECURVE_ERROR;
            }

            std::cout << "WaveCurveFactory: leaving Rarefaction. Continue as: " << type[future_curve - 1] << std::endl;

        } // if (current_curve == RAREFACTION_CURVE)
        else if (current_curve == COMPOSITE_CURVE){
            std::cout << "WaveCurveFactory: entering Composite." << std::endl;

//            double deltaxi = 1e-3;
            double deltaxi = 1e-3; //3e-4; // Was: 1e-3

            Curve cmpcurve, new_rarcurve;
            std::vector<int> index_explicit_bifurcation_transition;

            RealVector final_direction;

            int composite_stopped_because;

            int info_cmp = compositecurve->curve(g, f, b,
                                                 hwc.wavecurve[rarefaction_list.back()],
                                                 hwc.wavecurve[rarefaction_list.back()].curve[last_point_in_rarefaction.back()],
                                                 last_point_in_rarefaction.back(), 
                                                 odesolver, deltaxi,
                                                 linobj, linear_function,
                                                 future_composite_type, // COMPOSITE_BEGINS_AT_INFLECTION or COMPOSITE_AFTER_COMPOSITE.
                                                 family, 
                                                 new_rarcurve,
                                                 cmpcurve,
                                                 final_direction,
                                                 composite_stopped_because,
                                                 edge);

            is_first = false;


            std::cout << "WaveCurveFactory, composite completed. Info = " << info_cmp << ", final_direction = " << final_direction << std::endl;
            cmpcurve.final_direction = final_direction;

//            hwc.wavecurve.push_back(cmpcurve);

            if (info_cmp == COMPOSITE_OK){
                cmpcurve.back_curve_index = rarefaction_list.back();

                cmpcurve.back_curve_pointer.resize(cmpcurve.curve.size());
                for (int i = 0; i < cmpcurve.curve.size(); i++) cmpcurve.back_curve_pointer[i] = rarefaction_list.back();

                std::cout << "cmpcurve.back_curve_index = " << cmpcurve.back_curve_index << std::endl;
                std::cout << "composite_stopped_because = " << composite_stopped_because << std::endl;

                
                hwc.wavecurve.push_back(cmpcurve);

                future_curve_initial_point     = cmpcurve.last_point;
                future_curve_initial_direction = cmpcurve.final_direction;

                if (composite_stopped_because == COMPOSITE_REACHED_DOUBLE_CONTACT){
                    future_curve = RAREFACTION_CURVE;

                    std::cout << "WaveCurveFactory. After CompositeCurve reached a double contact." << std::endl;

                    last_point_in_rarefaction.back() = cmpcurve.back_pointer.back() - 1; // CHECK THIS!!!! Could be +/- 1 or something!
                    std::cout << "    last_point_in_rarefaction.back() = " << last_point_in_rarefaction.back() << std::endl;

                    if (last_point_in_rarefaction.back() < 0){ 
                        std::cout << "    Error! Memory positions cannot be negative! Aborting now." << std::endl;
                        return WAVECURVE_ERROR;
                    }
                }
                else if (composite_stopped_because == COMPOSITE_COMPLETED){
                    rarefaction_list.pop_back();
                    last_point_in_rarefaction.pop_back();

                    if (rarefaction_list.size() > 0){
                        future_curve = COMPOSITE_CURVE;
                        future_composite_type = COMPOSITE_AFTER_COMPOSITE;

                        compositecurve->correct_last_point(odesolver, deltaxi, hwc);
                    }
                    else { // Initial point reached!
                        future_curve = SHOCK_CURVE;
                    }
                }
                else if (composite_stopped_because == COMPOSITE_REACHED_BOUNDARY){
                    wavecurve_stopped_because = WAVECURVE_REACHED_BOUNDARY;

                    std::cout << "WaveCurveFactory: leaving Composite (reached boundary)." << std::endl;

                    return WAVECURVE_OK;
                }
                else if (composite_stopped_because == COMPOSITE_REACHED_LINE){
                    wavecurve_stopped_because = WAVECURVE_REACHED_LINE;

                    std::cout << "WaveCurveFactory: leaving Composite (reached line)." << std::endl;

                    return WAVECURVE_OK;
                }
            }
            else { // info_cmp == COMPOSITE_ERROR
                return  WAVECURVE_ERROR;
            }

            std::cout << "WaveCurveFactory: leaving Composite. Continue as: " << type[future_curve - 1] << std::endl;
        } // if (current_curve == COMPOSITE_CURVE)
        else if (current_curve == SHOCK_CURVE){
            std::cout << "WaveCurveFactory: entering Shock. direction = " << current_curve_initial_direction << std::endl;

            Curve shkcurve; 
            std::vector<int> stop_right_index;
            std::vector<int> stop_right_family;
            std::vector<int> stop_reference_index;    
            std::vector<int> stop_reference_family;  
        
            int shock_stopped_because;

            std::cout << "WaveCurveFactory: really entering Shock." << std::endl;

            // 
            int shck_info = shockcurve->curve_engine(hwc.reference_point, current_curve_initial_point, current_curve_initial_direction, family, 
                                            SHOCKCURVE_SHOCK_CURVE, 
                                            DONT_CHECK_EQUALITY_AT_LEFT /*SHOCK_SIGMA_EQUALS_LAMBDA_OF_FAMILY_AT_LEFT*/,
                                            SHOCK_SIGMA_EQUALS_LAMBDA_OF_FAMILY_AT_RIGHT,
                                            USE_ALL_FAMILIES, // TODO: Only for the time being! The correct one is USE_ALL_FAMILIES.
                                            STOP_AFTER_TRANSITION /*int after_transition*/,
                                            linobj, linear_function,
                                            shkcurve, 
                                            stop_right_index,
                                            stop_right_family,
                                            stop_reference_index,
                                            stop_reference_family,  
                                            shock_stopped_because,
                                            edge);

            // Kludge to solve the fact that the shockspeed at the reference point is being returned as NaN, 
            // which messes up with the Riemann Profile.
            //
            if (is_first){
                shkcurve.speed[0] = hwc.reference_point.e[family].r;
//                shkcurve.speed[0] = shkcurve.speed.back();
            }

            is_first = false;

            std::cout << "WaveCurveFactory. shck_info = " << shck_info << ", shock_stopped_because = " << shock_stopped_because << std::endl;

            std::cout << "Speed at first shockpoint = " << shkcurve.speed[0] << std::endl;


            shkcurve.back_curve_index = hwc.wavecurve.size() - 1;
            hwc.wavecurve.push_back(shkcurve);

            if (shck_info == SHOCKCURVE_OK){
                future_curve_initial_point = shkcurve.last_point;
                future_curve_initial_direction = shkcurve.final_direction;

                std::cout << "WaveCurveFactory, will try to store what the Shock produced." << std::endl;

                // TODO: This works sometimes, sometimes it doesn't. Why?
//                hwc.wavecurve.push_back(shkcurve);
//                hwc.add(shkcurve);

                std::cout << "WaveCurveFactory, successfully stored what the Shock produced." << std::endl;

                // The action depends whether forwards- or backwards-wavecurves are being constructed.

                if (shock_stopped_because == SHOCK_LEFT_CHARACTERISTIC_AT_FAMILY){
                    // Nothing to do in this case. The wavecurve stops.
                    // If other segments are desired continue the construction of an inadmissible Hugoniot.
                    // However, Panters suspects that this case does not arise if starting from a Lax shock.
                    // TODO: The statement above may be wrong. A shock may be invalid but may be followed by a valid one.


                    // TODO: In the future continue calculating the Hugoniot curve that represents an inadmissible shock.
                    //       Check for transitions that make it admissible again.

                    return WAVECURVE_OK;
                }
                else if (shock_stopped_because == SHOCK_RIGHT_CHARACTERISTIC_AT_FAMILY){
//                    std::vector<eigenpair> e;
//                    int n = shkcurve.curve.back();
//                    JetMatrix F_J(n), G_J(n);
//                    f->jet(shkcurve.curve.back(), F_J, 1);
//                    g->jet(shkcurve.curve.back(), G_J, 1);

//                    Eigen::eig(n, F_J.Jacobian().data(), G_J.Jacobian().data(), e);

//                    for (int i = 0; i < e.size(); i++) {
//                        std::cout << "Family: " << i << std::endl;
//                        std::cout << "    lambda = " << e[i].r << std::endl;
//                        std::cout << "    r = " << RealVector(n, e[i].vrr.data()) << std::endl << std::endl;
//                    }

//                    TestTools::pause("Check console.");

//                    std::cout << "Transitions: " << stop_current_family.size() << std::endl;
//                    for (int i = 0; i < stop_current_family.size(); i++) std::cout << "    " << stop_current_family[i] << std::endl;
//                    TestTools::pause("Check console.");

                    // TODO: In the future continue calculating the Hugoniot curve that represents an inadmissible shock.
                    //       Check for transitions that make it admissible again.

                    if (stop_right_family.size() > 0) {
//                        if (stop_right_family.back() == family) future_curve = RAREFACTION_CURVE;
//                        else return WAVECURVE_OK;

                        family = stop_right_family.back();
                        future_curve = RAREFACTION_CURVE;
                    }
                    else {
                        return WAVECURVE_OK;
                    }
                }
//                else if (shock_stopped_because == SHOCK_LEFT_CHARACTERISTIC_AT_CONTIGUOUS_FAMILY){
//                       
//                }
//                else if (shock_stopped_because == SHOCK_RIGHT_CHARACTERISTIC_AT_CONTIGUOUS_FAMILY){
//                else if (shock_stopped_because == SHOCK_RIGHT_CHARACTERISTIC_AT_UPPER_FAMILY){
//                    // Panters: The classical wavecurve must stop here. 
//                    // There is not a constant state between the wavecurves of the current and the upper family.

//                    return WAVECURVE_OK;                    
//                }
//                else if (shock_stopped_because == SHOCK_LEFT_CHARACTERISTIC_AT_LOWER_FAMILY){
//                    // If the last curve of the lower family's wavecurve was a rarefaction 
//                    // there is not a constant state between the wavecurves of both the current and the lower family.

//                    return WAVECURVE_OK;
//                }
//                else if (shock_stopped_because == SHOCK_LEFT_CHARACTERISTIC_AT_CONTIGUOUS_FAMILY){
//                    // TODO: Study this matter with Vitor. Build an auxiliary inadmissible curve.
//                    //
////                    future_curve = SHOCK_CURVE; // Panters: There is no certainty that the next curve is a shock.

//                    return WAVECURVE_OK;
//                }
                else if (shock_stopped_because == SHOCK_COMPLEX_EIGENVALUE_AT_FAMILY){
                    // TODO: Study this matter with Vitor. Build an auxiliary inadmissible curve.
                    //

                    // TODO: In the future continue calculating the Hugoniot curve that represents an inadmissible shock.
                    //       Check for transitions that make it admissible again.

                    future_curve = SHOCK_CURVE;

                    return WAVECURVE_OK;
                }
                else if (shock_stopped_because == SHOCK_REACHED_BOUNDARY){
                    wavecurve_stopped_because = WAVECURVE_REACHED_BOUNDARY;

                    

                    return WAVECURVE_OK;
                }
                else if (shock_stopped_because == SHOCK_REACHED_LINE){
                    wavecurve_stopped_because = WAVECURVE_REACHED_LINE;

                    std::cout << "WaveCurveFactory, shock reached line." << std::endl;

                    return WAVECURVE_OK;
                }
                else {

                    return WAVECURVE_OK;
                }
            }
            else {
                return WAVECURVE_ERROR;
            }

            std::cout << "WaveCurveFactory: leaving Shock. Continue as: " << type[future_curve - 1] << std::endl;

        } // if (current_curve == SHOCK_CURVE)
    }
}

// If the wavecurve starts at a boundary then it must call only ONE half_wavecurve.
//
// REMEMBER THAT THE RAREFACTION WILL BE IMPLEMENTED AS DIFFERENT CLASSES, AS HUGONIOT IS NOW.
// WHEN THE TIME COMES, THE RAREFACTION WILL BE PASSED AS A POINTER ALSO.
//
// TODO: Create the rarefaction, shock and composite externally, pass them here, use them as pointers in all the methods of WaveCurve.
//
int WaveCurveFactory::wavecurve(int type, const RealVector &initial_point, int family, int increase, HugoniotContinuation *h, 
                                void *linobj, double (*linear_function)(void *o, const RealVector &p),
                                WaveCurve &hwc, 
                                int &wavecurve_stopped_because, int &edge){

    // Initialize.
    //
    RealVector initial_direction;
    double dd;

    int info_initialize = rarefactioncurve->initialize(initial_point, family, increase, initial_direction, dd);

    if (info_initialize == RAREFACTION_INIT_ERROR) return WAVECURVEFACTORY_INIT_ERROR;

    ReferencePoint ref(initial_point, f, g, 0);

    hwc.family          = family;
    hwc.increase        = increase;
    hwc.reference_point = ref;

    // Proceed.
    //
    hugoniot = h;

    std::cout << "initial_direction = " << initial_direction << ", fam. = " << family << ", inc. = " << increase << std::endl;


    Liu_half_wavecurve(ref, initial_point, family, increase, RAREFACTION_CURVE,  initial_direction, linobj, linear_function, hwc, wavecurve_stopped_because, edge);

    hwc.beginnig_of_second_half = hwc.wavecurve.size();

    Liu_half_wavecurve(ref, initial_point, family, increase, SHOCK_CURVE,       -initial_direction, linobj, linear_function, hwc, wavecurve_stopped_because, edge);

    for (int i = 0; i < hwc.wavecurve.size(); i++) std::cout << "Curve\'s size = " << hwc.wavecurve[i].curve.size() << std::endl;

    add_arclength(0, hwc.beginnig_of_second_half - 1, 1.0, hwc);
    add_arclength(hwc.beginnig_of_second_half, hwc.wavecurve.size() - 1, -1.0, hwc);

    return WAVECURVE_OK;
}

void WaveCurveFactory::add_arclength(int begin, int end, double factor, WaveCurve &hwc){
    double distance = 0.0;
    RealVector prev_point = hwc.wavecurve[begin].curve.front();

    for (int i = begin; i <= end; i++){
        hwc.wavecurve[i].xi.resize(hwc.wavecurve[i].curve.size());
        for (int j = 0; j < hwc.wavecurve[i].curve.size(); j++){
            distance += norm(hwc.wavecurve[i].curve[j] - prev_point);
            hwc.wavecurve[i].xi[j] = distance*factor;

            prev_point = hwc.wavecurve[i].curve[j];
        }
    }

    return;
}

int WaveCurveFactory::wavecurve_from_boundary(const RealVector &initial_point, int s, int family, int increase, HugoniotContinuation *h, WaveCurve &hwc, int &wavecurve_stopped_because, int &edge){
    // Points to the interior of the domain from side s.
    //
    RealVector to_interior = b->side_transverse_interior(initial_point, s);

    // Initialize.
    //
    RealVector initial_direction;
    double dd;

    int info_initialize = rarefactioncurve->initialize(initial_point, family, increase, initial_direction, dd);

    if (info_initialize == RAREFACTION_INIT_ERROR) return WAVECURVEFACTORY_INIT_ERROR;

    ReferencePoint ref(initial_point, f, g, 0);

    hwc.family          = family;
    hwc.increase        = increase;
    hwc.reference_point = ref;

    // Proceed.
    //
    hugoniot = h;

    // PLEASE NOTE:
    // The case when initial_direction*to_interior == 0 is treated here by default. We do not know how to proceed correctly
    // in this case.
    //
    if (initial_direction*to_interior > 0.0){
        Liu_half_wavecurve(ref, initial_point, family, increase, RAREFACTION_CURVE,  initial_direction, hwc, wavecurve_stopped_because, edge);
    }
    else {
        Liu_half_wavecurve(ref, initial_point, family, increase, SHOCK_CURVE,       -initial_direction, hwc, wavecurve_stopped_because, edge);    
    }

    return WAVECURVE_OK;
}

int WaveCurveFactory::wavecurve_from_inflection(const std::vector<RealVector> &inflection_curve, const RealVector &p, int family, int increase, HugoniotContinuation *h, WaveCurve &hwc, int &wavecurve_stopped_because, int &edge){
    // Proceed.
    //
    hugoniot = h;

    family_for_directional_derivative = family;

    RealVector closest_point;
    Utilities::pick_point_from_segmented_curve(inflection_curve, p, closest_point);

    reference_for_directional_derivative = closest_point - p;
    normalize(reference_for_directional_derivative);

    RealVector point_on_level_curve;
    int info_find_point_on_level_curve = Utilities::find_point_on_level_curve((void*)this, &rarefaction_directional_derivative, p, closest_point, point_on_level_curve);

    // Compute the eigenvectors in point_on_level_curve.
    //
    std::vector<eigenpair> e;
    Eigen::fill_eigenpairs(f, g, point_on_level_curve, e);

    double deltaxi = 1e-3;
    int n = p.size();

    std::vector<RealVector> direction;
    direction.push_back(RealVector(n, e[family].vrr.data()));
    direction.push_back(-direction[0]);

    ReferencePoint ref(point_on_level_curve, f, g, 0);

    hwc.family          = family;
    hwc.increase        = increase;
    hwc.reference_point = ref;

    for (int i = 0; i < direction.size(); i++){
        std::vector<double> lambda;
        Eigen::fill_eigenvalues(f, g, point_on_level_curve + deltaxi*direction[i], lambda);

        int start_as;

        std::cout << "i = " << i << ", lambda = " << lambda[family] << ", e[family].r = " << e[family].r << std::endl;


        if (
            (lambda[family] > e[family].r && increase == SPEED_INCREASE) ||
            (lambda[family] < e[family].r && increase == SPEED_DECREASE)
           ) {
            start_as = RAREFACTION_CURVE;
        }
        else {
            start_as = SHOCK_CURVE;
        }

        std::cout << "start as: " << start_as << std::endl; 
        Liu_half_wavecurve(ref, point_on_level_curve + deltaxi*direction[i], family, increase, start_as, direction[i], hwc, wavecurve_stopped_because, edge);
    }

    return WAVECURVE_OK;
}

double WaveCurveFactory::rarefaction_directional_derivative(void *obj, const RealVector &p){
    WaveCurveFactory *wavecurve = (WaveCurveFactory*)obj;

    RarefactionCurve *rarefactioncurve = wavecurve->rarefactioncurve;
    int family = wavecurve->family_for_directional_derivative; 
    RealVector reference_for_directional_derivative = wavecurve->reference_for_directional_derivative;

    return rarefactioncurve->directional_derivative(p, family, reference_for_directional_derivative);
}

int WaveCurveFactory::intersection(const WaveCurve &c1, const WaveCurve &c2, const RealVector &pmin, const RealVector &pmax, 
                                   RealVector &p, int &subc1, int &subc1_point, int &subc2, int &subc2_point){

    // This will be removed in the near future. The wavecurves will have an accompanying hyperoctree. Thus,
    // there will be no need to build said hyperoctree here. Alternatively, the domain will not be [0, 0]-[1, 1]
    // and therefor the lines below will also be wrong.
    PointND hpmin(2), hpmax(2);
    hpmin(0) = hpmin(1) = 0.0;
    hpmax(0) = hpmax(1) = 1.0;
    BoxND b(hpmin, hpmax);

    HyperOctree<WaveCurveSegment> h1(b), h2(b);

    // Fill quadtrees
    for (int i = 0; i < c1.wavecurve.size(); i++){
        if (c1.wavecurve[i].curve.size() < 2) continue;
        for (int j = 0; j < c1.wavecurve[i].curve.size() - 1; j++){
            //std::cout << "1, (" <<  c1[i].curve[j] << ")-(" << c1[i].curve[j + 1] << ")" << std::endl;
            WaveCurveSegment *wcs = new WaveCurveSegment(c1.wavecurve[i].curve[j], c1.wavecurve[i].curve[j + 1], i, j);
            h1.add(wcs);
        }
    }

    for (int i = 0; i < c2.wavecurve.size(); i++){
        if (c2.wavecurve[i].curve.size() < 2) continue;
        for (int j = 0; j < c2.wavecurve[i].curve.size() - 1; j++){
            //std::cout << "2, (" <<  c2[i].curve[j] << ")-(" << c2[i].curve[j + 1] << ")" << std::endl;
            WaveCurveSegment *wcs = new WaveCurveSegment(c2.wavecurve[i].curve[j], c2.wavecurve[i].curve[j + 1], i, j);
            h2.add(wcs);
        }
    }

    // Query the quadtrees to find the segments contained in the given box.
    PointND ppmin(2), ppmax(2);            // Add lambda as the third component of each point of the composite curve

    for (int i = 0; i < 2; i++){
        ppmin(i) = pmin.component(i);
        ppmax(i) = pmax.component(i);
    }

    BoxND box(ppmin, ppmax);

    std::vector<WaveCurveSegment*> wcs1, wcs2;
    h1.within_box(box, wcs1);
    h2.within_box(box, wcs2);

    // Find the first intersection point:
    //
    bool found = false;

    for (int i = 0; i < wcs1.size(); i++){
        for (int j = 0; j < wcs2.size(); j++){
            RealVector r(2);

            if (segment_intersection(wcs1[i]->rv[0].components(), wcs1[i]->rv[1].components(), 
                                     wcs2[j]->rv[0].components(), wcs2[j]->rv[1].components(), 
                                     r.components())
               ) {
                p = r;
                subc1 = wcs1[i]->curve_position; subc1_point = wcs1[i]->segment_position;
                subc2 = wcs2[j]->curve_position; subc2_point = wcs2[j]->segment_position;

                found = true;
                continue;
            }
        }
    }

    for (int i = 0; i < wcs1.size(); i++) delete wcs1[i];
    for (int i = 0; i < wcs2.size(); i++) delete wcs2[i];

    if (found) return WAVE_CURVE_INTERSECTION_FOUND;
    else       return WAVE_CURVE_INTERSECTION_NOT_FOUND;
}

int WaveCurveFactory::wavecurve_from_wavecurve(const WaveCurve &c, const RealVector &p, HugoniotContinuation *h, WaveCurve &hwc, int &wavecurve_stopped_because, int &edge){

    int curve_index, segment_index_in_curve;
    RealVector closest_point;
    double speed;

    Utilities::pick_point_from_wavecurve(c, p, curve_index, segment_index_in_curve, closest_point, speed);

    int family = (c.increase == SPEED_INCREASE) ? c.family + 1 : c.family - 1;

    return wavecurve(WAVECURVEFACTORY_GENERIC_POINT, closest_point, family, c.increase, h, hwc, wavecurve_stopped_because, edge);
}

//int WaveCurveFactory::wavecurve_from_wavecurve(const WaveCurve &c, const RealVector &p, HugoniotContinuation *h, WaveCurve &hwc, int &wavecurve_stopped_because, int &edge){

//    int curve_index, segment_index_in_curve;
//    RealVector closest_point;
//    double speed;

//    Utilities::pick_point_from_wavecurve(c, p, curve_index, segment_index_in_curve, closest_point, speed);

////    int family = (c.increase == SPEED_INCREASE) ? c.family + 1 : c.family - 1;

//    int increase = (c.increase == SPEED_INCREASE) ? SPEED_DECREASE : SPEED_INCREASE;

//    return wavecurve(closest_point, c.family, increase, h, hwc, wavecurve_stopped_because, edge);
//}

//int WaveCurveFactory::wavecurve_from_right_state(const WaveCurve &c, const RealVector &p, HugoniotContinuation *h, WaveCurve &hwc, int &wavecurve_stopped_because, int &edge){

//    int curve_index, segment_index_in_curve;
//    RealVector closest_point;
//    double speed;

//    Utilities::pick_point_from_wavecurve(c, p, curve_index, segment_index_in_curve, closest_point, speed);

//    int family = (c.increase == SPEED_INCREASE) ? c.family + 1 : c.family - 1;

//    return wavecurve(closest_point, family, c.increase, h, hwc, wavecurve_stopped_because, edge);
//}


void WaveCurveFactory::R_regions(HugoniotContinuation *h, const WaveCurve &c, std::vector<WaveCurve> &curves){
    curves.clear();

    int increase = c.increase;
    int family = (c.increase == SPEED_INCREASE) ? c.family + 1: c.family - 1;

    for (int i = 1; i < c.wavecurve.size(); i++){
        WaveCurve w;

        int wavecurve_stopped_because, edge;

        wavecurve(WAVECURVEFACTORY_GENERIC_POINT, c.wavecurve[i].curve[0], family, increase, h, w, 
                  wavecurve_stopped_because, edge);

        curves.push_back(w);
    }

    return;
}


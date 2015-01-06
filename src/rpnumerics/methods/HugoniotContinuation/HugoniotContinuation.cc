#include "HugoniotContinuation.h"


// Trivial initialization. It is assumed that the reference point and the initial point
// are one and the same (in).
//
// TODO: Accumulations should be able to declare if they are trivial or not.
//       If a trivial accumulation is used, see below to modify some parts accordingly.
//
int HugoniotContinuation::find_initial_direction(const RealVector &in, const RealVector &hint_direction, int fam, RealVector &initial_direction){
    int n = in.size();

    DoubleMatrix A(n, n), B(n, n);
    f->fill_with_jet(n, in.components(), 1, 0, A.data(), 0);
    g->fill_with_jet(n, in.components(), 1, 0, B.data(), 0);
    
    RealVector r;

    // TODO: If the accumulation is trivial this method should
    //       be replaced by the one that only takes one matrix as input.
    //
    int info;
    info = Eigen::eig(n, A.data(), B.data(), fam, r); // TODO: Verify if lambda is complex, return error. Make similar for eigenvalue. Error not fatal in that case.
    
    //if (g != 0) info = Eigen::eig(n, A.data(), B.data(), fam, r);
    //else        info = Eigen::eig(n, A.data(), fam, r);

    // TODO: Pensar como corrigir de modo a tornar hyperplane paralelo aos autovetores distintos do de fam. (Dan)

    if (info == COMPLEX_EIGENVALUE) return HUGONIOTCONTINUATION_INITIALIZE_ERROR; 

    if (r*hint_direction > 0.0) initial_direction = r;
    else                        initial_direction = -r;

    return HUGONIOTCONTINUATION_INITIALIZED_OK;
}

// For the case when the reference point is not point where the continuation is found.
// Also, when the wavecurve is computed, this method will be used to compute a shock.
// Therefore, this method must be public.
//
int HugoniotContinuation::find_continuation_direction(const RealVector &Hugoniot_point, const RealVector &hint, RealVector &Hugoniot_direction){

    DoubleMatrix hyperplane;
    int info_fill_hyperplane = fill_hyperplane(Hugoniot_point, hyperplane);
    if (info_fill_hyperplane != HUGONIOTCONTINUATION_HYPERPLANE_OK) return info_fill_hyperplane;
    
    int info_fill_Hugoniot_direction = fill_Hugoniot_direction(hint, hyperplane, Hugoniot_direction);
    if (info_fill_Hugoniot_direction != HUGONIOTCONTINUATION_DIRECTION_OK) return info_fill_Hugoniot_direction;
    else return HUGONIOTCONTINUATION_DIRECTION_OK;
}

// Shockspeed. Not used right now.
// 
double HugoniotContinuation::shockspeed(const FluxFunction *fp, const AccumulationFunction *gp, const RealVector &Up,
                                        const FluxFunction *fm, const AccumulationFunction *gm, const RealVector &Um) const {

    int n = Up.size();

    RealVector Fp(n), Fm(n), Gp(n), Gm(n);

    fp->fill_with_jet(n, Up.components(), 0, Fp.components(), 0, 0);
    fm->fill_with_jet(n, Um.components(), 0, Fm.components(), 0, 0);

    gp->fill_with_jet(n, Up.components(), 0, Gp.components(), 0, 0);
    gm->fill_with_jet(n, Um.components(), 0, Gm.components(), 0, 0);

    double s = 0.0, d = 0.0;

    for (int i = 0; i < n; i++){
        double dd = Gp(i) - Gm(i);

        s += (Fp(i) - Fm(i))*dd;

        d += dd*dd;
    }

    return s/d; 
}

void HugoniotContinuation::jet_sigma(const RealVector &U, const RealVector &Hugoniot_direction, double &sigma, double &sigma_dot) const {
    int n = ref.point.size();

    RealVector F(n), G(n);
    DoubleMatrix JF(n, n), JG(n, n);

    f->fill_with_jet(n, U.components(), 1, F.components(), JF.data(), 0);
    g->fill_with_jet(n, U.components(), 1, G.components(), JG.data(), 0);

    RealVector diff_F = F - ref.F;
    RealVector diff_G = G - ref.G;

    double inv_diff_G_norm_squared = 1.0/(diff_G*diff_G);

    sigma = (diff_F*diff_G)*inv_diff_G_norm_squared;

    RealVector P = (JF - sigma*JG)*Hugoniot_direction;

    sigma_dot = (P*diff_G)*inv_diff_G_norm_squared;

    return;
}

// Sigma     = s(0)
// Sigma_dot = s(0, 0)
//
void HugoniotContinuation::sigma_jet(const RealVector &U, const RealVector *Hugoniot_direction, int degree, JetMatrix &s) const {
    s.resize(1);
    
    if (degree >= 0){
        int n = ref.point.size();

        RealVector F(n), G(n);
        DoubleMatrix JF(n, n), JG(n, n);

        f->fill_with_jet(n, U.components(), degree, F.components(), JF.data(), 0);
        g->fill_with_jet(n, U.components(), degree, G.components(), JG.data(), 0);

        RealVector diff_F = F - ref.F;
        RealVector diff_G = G - ref.G;

        double inv_diff_G_norm_squared = 1.0/(diff_G*diff_G);

        double sigma = (diff_F*diff_G)*inv_diff_G_norm_squared; 
        s.set(0, sigma);
        
        if (degree >= 1){
            RealVector P = (JF - sigma*JG)*(*Hugoniot_direction);

            double sigma_dot = (P*diff_G)*inv_diff_G_norm_squared;
            s.set(0, 0, sigma_dot);
        }
    }

    return;
}

double HugoniotContinuation::sigma(const RealVector &F, const RealVector &G) const {
    RealVector diff_F = F - ref.F;
    RealVector diff_G = G - ref.G;

    double inv_diff_G_norm_squared = 1.0/(diff_G*diff_G);
    
    return (diff_F*diff_G)*inv_diff_G_norm_squared;
}

double HugoniotContinuation::sigma(const RealVector &Fp, const RealVector &Gp, const RealVector &Fm, const RealVector &Gm) const {
    RealVector diff_F = Fp - Fm;
    RealVector diff_G = Gp - Gm;

    double inv_diff_G_norm_squared = 1.0/(diff_G*diff_G);
    
    return (diff_F*diff_G)*inv_diff_G_norm_squared;
}

int HugoniotContinuation::fill_hyperplane(const RealVector &origin, DoubleMatrix &hyperplane){
//    int n = origin.size();

//    RealVector F(n), G(n);
//    DoubleMatrix JF(n, n), JG(n, n);

//    f->fill_with_jet(n, origin.components(), 1, F.components(), JF.data(), 0);
//    g->fill_with_jet(n, origin.components(), 1, G.components(), JG.data(), 0);
//    jet_Hugoniot(origin, H, hyperplane);
//    jet_Hugoniot(F, JF, G, JG, H, hyperplane);

    JetMatrix Fjet, Gjet;
    f->jet(origin, Fjet, 1);
    g->jet(origin, Gjet, 1);

    RealVector H;
    jet_Hugoniot(Fjet.function(), Fjet.Jacobian(), Gjet.function(), Gjet.Jacobian(), H, hyperplane);

    // Normalize by cols.
    //
    for (int i = 0; i < hyperplane.cols(); i++){
        double nrm = 0.0;
        for (int j = 0; j < hyperplane.rows(); j++) nrm += hyperplane(j, i)*hyperplane(j, i);

        double inv_nrm = 1.0/sqrt(nrm);

        for (int j = 0; j < hyperplane.rows(); j++) hyperplane(j, i) *= inv_nrm;
    }

    return HUGONIOTCONTINUATION_HYPERPLANE_OK;
}

void HugoniotContinuation::Newton_system_for_Hugoniot(const RealVector &p, const DoubleMatrix &hyperplane, DoubleMatrix &Newton_matrix, RealVector &error){
    DoubleMatrix nablaH;

//    int n = p.size();

//    RealVector F(n), G(n);
//    DoubleMatrix JF(n, n), JG(n, n);

//    f->fill_with_jet(n, p.components(), 1, F.components(), JF.data(), 0);
//    g->fill_with_jet(n, p.components(), 1, G.components(), JG.data(), 0);

    JetMatrix Fjet, Gjet;
    f->jet(p, Fjet, 1);
    g->jet(p, Gjet, 1);
    
//    jet_Hugoniot(p, error, nablaH);
//    jet_Hugoniot(F, JF, G, JG, error, nablaH);
    jet_Hugoniot(Fjet.function(), Fjet.Jacobian(), Gjet.function(), Gjet.Jacobian(), error, nablaH);

    Newton_matrix = transpose(nablaH)*hyperplane;

    return;
}

int HugoniotContinuation::Newton_in_hyperplane(const RealVector &origin, const DoubleMatrix &hyperplane, RealVector &Hugoniot_intersection){

    RealVector correction(RealVector::zeroes(origin.size() - 1));
    
    int max_number_iterations = 10; // TODO: This parameter should receive a default value in the ctor.
    int iterations = 0;
    
    double min_norm = 1e-6; // TODO: This parameter should receive a default value in the ctor.
    
    double norm_correction = 0.0;
    
    RealVector hyperplane_point = origin;
    
    do {
        DoubleMatrix Newton_matrix;
        RealVector error;

        Newton_system_for_Hugoniot(hyperplane_point, hyperplane, Newton_matrix, error);

        int info_correction = solve(Newton_matrix, error, correction);
        
        if (info_correction == REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR) return REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR;
        
        // Update the point in the hyperplane.
        //
        hyperplane_point = hyperplane_point - hyperplane*correction;

        iterations++;
        
        norm_correction = norm(correction);
        
    } while (iterations < max_number_iterations && norm_correction > min_norm && norm_correction < 10.0); // TODO: 10.0 was fixed here, it must be fine-tuned, perhaps it should become a member.
    
    if (iterations >= max_number_iterations) return HUGONIOTCONTINUATION_NEWTON_ERROR;
    
    // Output.
    //
    Hugoniot_intersection = hyperplane_point;
    
    return HUGONIOTCONTINUATION_NEWTON_OK;
}

// TODO: Chamar rotina para evitar cruzar secondary bifurcation entre previous hyperplane origin e ne hyperplane origin.
// Corrigir a direcao do hyperplane, talvez. (Dan)
//
int HugoniotContinuation::Newton_step(const RealVector &previous_point, double &step_size, int &number_of_steps_with_unchanged_size,
                                      const RealVector &previous_direction,
                                      /*DoubleMatrix &hyperplane, */RealVector &Hugoniot_intersection){
    int step_size_iterations = 0;
    int max_number_step_size_iterations = 10;
    int info_newton;
    double cos_angle = 1.0;

    RealVector hyperplane_origin;

    step_size *= 2.0; 

    DoubleMatrix hyperplane;

    do {
        step_size *= 0.5;
        step_size_iterations++;

        hyperplane_origin = previous_point + step_size*previous_direction;

//        std::cout << "hyperplane_origin = " << hyperplane_origin << std::endl;

        // The hyperplane where the Hugoniot Locus point will be found via Newton.
        // (The columns are the basis' vectors).
        //
        int info_fill_hyperplane = fill_hyperplane(hyperplane_origin, hyperplane);

//        std::cout << "hyperplane = " << hyperplane << std::endl;
        
        if (info_fill_hyperplane == HUGONIOTCONTINUATION_HYPERPLANE_ERROR){
            return HUGONIOTCONTINUATION_HYPERPLANE_ERROR;
        } 
        
        // Find the intersection of the Hugoniot curve with the hyperplane.
        //
        info_newton = Newton_in_hyperplane(hyperplane_origin, hyperplane, Hugoniot_intersection);

        // Check that the angle between previous_direction and (Hugoniot_intersection - previous_point) does not exceed arccos(MAX_COS_ANGLE). 
        // If so, reduce shift.
        //
        RealVector new_candidate_direction = Hugoniot_intersection - previous_point;
        normalize(new_candidate_direction);

        cos_angle = new_candidate_direction*previous_direction;

    } while (step_size_iterations < max_number_step_size_iterations && (info_newton == HUGONIOTCONTINUATION_NEWTON_ERROR || cos_angle <= MAX_COS_ANGLE) );

    // If only one step was needed, increase number_of_steps_with_unchanged_size.
    // Hopefully step_size will increase also in the near future.
    //
    if (step_size_iterations == 1) number_of_steps_with_unchanged_size++;
    else                           number_of_steps_with_unchanged_size = 0;

    // TODO: Here, if convergence is not achieved, the shift should be
    //       reduced. When the moment comes this region will go
    //       inside a do-while loop.
    //
    if (info_newton != HUGONIOTCONTINUATION_NEWTON_OK){ // Was == HUGONIOTCONTINUATION_NEWTON_ERROR
        return HUGONIOTCONTINUATION_NEWTON_ERROR;
    } 
    else return HUGONIOTCONTINUATION_NEWTON_STEP_OK;
}

int HugoniotContinuation::Newton_step(const RealVector &previous_point, double &step_size, int &number_of_steps_with_unchanged_size,
                                      const RealVector &previous_direction,
                                      DoubleMatrix &hyperplane, RealVector &Hugoniot_intersection){
    int step_size_iterations = 0;
    int max_number_step_size_iterations = 10;
    int info_newton;
    double cos_angle = 1.0;

    RealVector hyperplane_origin;

    step_size *= 2.0; 

    do {
        step_size *= 0.5;
        step_size_iterations++;

        hyperplane_origin = previous_point + step_size*previous_direction;

        // The hyperplane where the Hugoniot Locus point will be found via Newton.
        // (The columns are the basis' vectors).
        //
        int info_fill_hyperplane = fill_hyperplane(hyperplane_origin, hyperplane);

        if (info_fill_hyperplane == HUGONIOTCONTINUATION_HYPERPLANE_ERROR){
            return HUGONIOTCONTINUATION_HYPERPLANE_ERROR;
        } 
        
        // Find the intersection of the Hugoniot curve with the hyperplane.
        //
        info_newton = Newton_in_hyperplane(hyperplane_origin, hyperplane, Hugoniot_intersection);

        // Check that the angle between previous_direction and (Hugoniot_intersection - previous_point) does not exceed arccos(MAX_COS_ANGLE). 
        // If so, reduce shift.
        //
        RealVector new_candidate_direction = Hugoniot_intersection - previous_point;
        normalize(new_candidate_direction);

        cos_angle = new_candidate_direction*previous_direction;

    } while (step_size_iterations < max_number_step_size_iterations && (info_newton == HUGONIOTCONTINUATION_NEWTON_ERROR || cos_angle <= MAX_COS_ANGLE) );

    // If only one step was needed, increase number_of_steps_with_unchanged_size.
    // Hopefully step_size will increase also in the near future.
    //
    if (step_size_iterations == 1) number_of_steps_with_unchanged_size++;
    else                           number_of_steps_with_unchanged_size = 0;

    // TODO: Here, if convergence is not achieved, the shift should be
    //       reduced. When the moment comes this region will go
    //       inside a do-while loop.
    //
    if (info_newton != HUGONIOTCONTINUATION_NEWTON_OK){ // Was == HUGONIOTCONTINUATION_NEWTON_ERROR
        return HUGONIOTCONTINUATION_NEWTON_ERROR;
    } 
    else return HUGONIOTCONTINUATION_NEWTON_STEP_OK;
}

int HugoniotContinuation::fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &hyperplane, RealVector &Hugoniot_direction){
    Hugoniot_direction = orthogonalize(previous_direction, hyperplane);
    
    // Normalize Hugoniot_direction and verify that is coherent with previous_direction.
    //
    double norm_Hugoniot_direction = norm(Hugoniot_direction);
    
    // TODO: The 1e-4 below should be variable and should come from outside somehow.
    //
    if (norm_Hugoniot_direction < 1e-4) return HUGONIOTCONTINUATION_DIRECTION_ERROR;
    else {
        double inv = 1.0/norm_Hugoniot_direction;
        Hugoniot_direction = Hugoniot_direction*inv;
    
        if (Hugoniot_direction*previous_direction < 0.0) Hugoniot_direction = -Hugoniot_direction;

        return HUGONIOTCONTINUATION_DIRECTION_OK;
    }
}

// TODO: When the curve is initialized the initial point should be added to the curve.
//
int HugoniotContinuation::curve_engine(const RealVector &in, const RealVector &initial_direction, 
                                       RealVector &final_direction, std::vector<RealVector> &shockcurve, int &edge){

    double step_size;

//    // This is the distance between two consecutive hyperplanes.
//    //
//    double default_step_size = 0.00005; /*0.005*/ // The order of magnitude of f and g must be known, so as to declare a correct value for default_shift.
//    double max_default_step_size = 100.0*default_step_size;

    RealVector previous_point = in;
    RealVector previous_direction = initial_direction;

    step_size = default_step_size_;
    int number_of_steps_with_unchanged_size = 0;
    int step_size_increased = 0;

    while (true){
        RealVector Hugoniot_intersection;
        RealVector Hugoniot_direction;

        // TODO: Decide which version of the following methods is to survive. Right now (16-Aug-2013) the version that doesn't write on hyperplane
        //       works fine.
        DoubleMatrix hyperplane;
//        int info_Newton_step = Newton_step(previous_point, step_size, number_of_steps_with_unchanged_size, previous_direction, hyperplane, Hugoniot_intersection);

        int info_Newton_step = Newton_step(previous_point, step_size, number_of_steps_with_unchanged_size, previous_direction, /*hyperplane, */Hugoniot_intersection);

        // Check if size_steps remains the same after a few steps.
        // If so, increase step_size;
        //
        if (number_of_steps_with_unchanged_size == 5){
            number_of_steps_with_unchanged_size = 0;
            
            // TODO: See if something is missing here.
            step_size_increased++;

            if (step_size_increased < 16){ // Was: 8
                step_size = min(step_size*SQRT_TWO, max_default_step_size_);
            }
        }
        
        if (info_Newton_step != HUGONIOTCONTINUATION_NEWTON_STEP_OK){
            return info_Newton_step;
        } 
        
        // Find the hyperplane.
        //
        fill_hyperplane(Hugoniot_intersection, hyperplane);
        
        // Find the direction of the Hugoniot curve.
        //
        int info_fill_Hugoniot_direction = fill_Hugoniot_direction(previous_direction, hyperplane, Hugoniot_direction);

        if (info_fill_Hugoniot_direction == HUGONIOTCONTINUATION_DIRECTION_ERROR){
            Hugoniot_direction = previous_direction;
        }
        
        // Something's missing here...Newton_step, leaving.

        // If the bifurcation space was crossed, add the point where the space was intersected (found by linear interpolation)
        // to the curve.

        previous_point     = Hugoniot_intersection; std::cout << Hugoniot_intersection << std::endl;
        previous_direction = Hugoniot_direction;
        
        // Verify if the new point lies within the domain or not.
        // TODO: Why is it that this is not working correctly for the Liquid?
        //
        if (shockcurve.size() > 1){
            RealVector r;
            int info_intersect = b->intersection(shockcurve[shockcurve.size() - 1], Hugoniot_intersection, r, edge);

            // Both points are inside: carry on.
            if (info_intersect == 1)       shockcurve.push_back(Hugoniot_intersection);
            // Both outside (this really should not happen).
            else if (info_intersect == -1) return HUGONIOTCONTINUATION_CURVE_ERROR;
            // New point outside: the curve reached the domain's boundary.
            else { 
                shockcurve.push_back(r);

                return HUGONIOTCONTINUATION_CURVE_OK;
            }
        }
        else shockcurve.push_back(Hugoniot_intersection);
    }

    return HUGONIOTCONTINUATION_CURVE_OK;
}

// This method should be called from within a curve constructor.
// By itself it only advances one step.
// 
// The following parameters must be meaningfully filled before entering here and analised after leaving:
//
//     step_size_increased, step_size, number_of_steps_with_unchanged_size.
//
// TODO: Pass also F, G, that they need not be computed again here.
//
int HugoniotContinuation::curve_point(const RealVector &previous_point, double previous_sigma_between_points, 
                                      const RealVector &direction, 
                                      int &step_size_increased,
                                      double &step_size, int &number_of_steps_with_unchanged_size, 
                                      RealVector &Hugoniot_intersection, double &sigma_between_points,
                                      RealVector &Hugoniot_direction){
                                      
//    std::cout << "HugoniotContinuation. Entering curve_point, direction = " << direction << std::endl;

//    // TODO: These parameters must belong to the object, and SHOULD be made public. THEY ENTER, CANT STAY HERE
//    // This is the distance between two consecutive hyperplanes.
//    //
//    double default_step_size = 0.00005; // The order of magnitude of f and g must be known, so as to declare a correct value for default_shift.
//    double max_default_step_size = 100.0*default_step_size;

    int info_Newton_step = Newton_step(previous_point, step_size, number_of_steps_with_unchanged_size, direction, Hugoniot_intersection);

//    std::cout << "Hugoniot_intersection = " << Hugoniot_intersection << std::endl;

    sigma_between_points = shockspeed(f, g, previous_point, f, g, Hugoniot_intersection);
    if (std::abs(sigma_between_points - previous_sigma_between_points) < .1*(std::abs(sigma_between_points) + std::abs(previous_sigma_between_points))){

        // Check if size_steps remains the same after a few steps.
        // If so, increase step_size;
        //
        if (number_of_steps_with_unchanged_size == 5){
            number_of_steps_with_unchanged_size = 0;
            
            // TODO: See if something is missing here.
            step_size_increased++;

            if (step_size_increased < 16){ // Was: 8
                step_size = min(step_size*SQRT_TWO, max_default_step_size_);
            }
        }
    }
    else {
        // TODO: The step_size should be decreased here!
    }

        
    if (info_Newton_step != HUGONIOTCONTINUATION_NEWTON_STEP_OK){
        return info_Newton_step;
    } 

    // Find the hyperplane.
    //
    DoubleMatrix hyperplane;
    fill_hyperplane(Hugoniot_intersection, hyperplane);
        
    // Find the direction of the Hugoniot curve.
    //
    int info_fill_Hugoniot_direction = fill_Hugoniot_direction(direction, hyperplane, Hugoniot_direction);

    if (info_fill_Hugoniot_direction == HUGONIOTCONTINUATION_DIRECTION_ERROR){
//        Hugoniot_direction = previous_direction;
    }
        
    // Something's missing here...Newton_step, leaving.

    // If the bifurcation space was crossed, add the point where the space was intersected (found by linear interpolation)
    // to the curve.

//    previous_point     = Hugoniot_intersection;
//    previous_direction = Hugoniot_direction;

    return HUGONIOTCONTINUATION_NEWTON_STEP_OK; // Find a better return code.
}

// TODO: Subphysics may be included somehow in the reference point.
void HugoniotContinuation::set_reference_point(const ReferencePoint &r){
    ref = r;
    reference_point=r;
    return;
}

// ATTENTION!!!
//
// Theta_index depends on how the SubPhysics was coded!!!
// So far Theta_index = 1 (the second variable), but it can be any other, since there is no
// standard or agreed form to code SubPhysics' FluxFunction & AccumulationFunction.
//
// This method CANNOT be called before a reference point is set!
//
int HugoniotContinuation::set_bifurcation_space_coordinate(int Theta_index){
    if (Theta_index < 0) there_is_a_bifurcation_space = false; // This serves as a reset.
    else {
        if (Theta_index != 1) return HUGONIOTCONTINUATION_SET_BIFURCATION_PLANE_ERROR;

        bifurcation_space_coordinate = ref.point(Theta_index); // Get Theta.
        bifurcation_space_coordinate_index = Theta_index;
    
        // Not used right now:

//        int n = ref.point.size();
//        bifurcation_space_basis.resize(n, n - 1);
//        bifurcation_space_basis(0, 0) = 1.0;
//        bifurcation_space_basis(1, 0) = 0.0;
//        bifurcation_space_basis(2, 0) = 0.0;
//    
//        bifurcation_space_basis(0, 1) = 0.0;
//        bifurcation_space_basis(1, 1) = 0.0;
//        bifurcation_space_basis(2, 1) = 1.0;
    
        there_is_a_bifurcation_space = true;
    }    
    
    return HUGONIOTCONTINUATION_SET_BIFURCATION_PLANE_OK;
}

bool HugoniotContinuation::reference_point_near_coincidence(){

    if (ref.e.size() == 0) return false;
    else return fabs(ref.e[0].r - ref.e[1].r)/(1.0 + fabs(ref.e[0].r) + fabs(ref.e[1].r)) < 2e-3;
}

HugoniotContinuation::HugoniotContinuation(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb) : HugoniotCurve(ff, gg, bb){
    f = ff;
    g = gg;
    b = bb;

//    there_is_a_bifurcation_space = false;

    default_step_size_ = 1.0e-4 /*0.00005*/;
    max_default_step_size_ = 100.0*default_step_size_;

    info_ = std::string("HugoniotContinuation");
    method_ = EXPLICIT_HUGONIOT;
}

HugoniotContinuation::~HugoniotContinuation(){
}

// The initial point is not added yet.
//
int HugoniotContinuation::curve(std::vector< std::vector<RealVector> > &curves){
    std::vector<int> edges;
    cout<<"Algo"<< reference_point.point<<"                 "<<ref.point<< endl;
    int info = all_curves(curves, edges);

    return info;
}

int HugoniotContinuation::disconnected_curve(const RealVector &point_on_edge, int entry_edge, std::vector<RealVector> &curve, int &edge){
    return HUGONIOTCONTINUATION_CURVE_OK;
}

int HugoniotContinuation::all_curves(std::vector< std::vector<RealVector> > &curves, std::vector<int> &edges){
    cout <<"Em all"<<endl;
    int info = connected_curve(curves, edges);
    cout <<"ddfasdfasd"<<endl;
    std::vector<std::vector<RealVector> > disconnected_temp;
    std::vector<int> disconnected_edges;

    for (int i = 0; i < curves.size(); i++){
        std::vector<RealVector> temp;
        int edge;

        disconnected_curve(curves[i].back(), edges[i], temp, edge);

        if (temp.size() > 0){
            disconnected_temp.push_back(temp);
            disconnected_edges.push_back(edge);
        }
    }

    for (int i = 0; i < disconnected_temp.size(); i++){
        curves.push_back(disconnected_temp[i]);
        edges.push_back(disconnected_edges[i]);
    }

    return info;
}

int HugoniotContinuation::connected_curve(std::vector< std::vector<RealVector> > &curves, std::vector<int> &edges){
    curves.clear();
    edges.clear();

    int n = ref.point.size();

    // If the initial point lies near the coincidence curve, abort.
    // TODO: These lines below need to be improved.
    //
    

    if (reference_point_near_coincidence()){ // TODO: This value must be fine-tuned.

//        initial_step_size = ...;
        return HUGONIOTCONTINUATION_NEAR_COINCIDENCE_CURVE;
    }

    for (int family = 0; family < ref.e.size(); family++){
        // Find the eigenvector of this family
        RealVector r(n);
        for (int i = 0; i < n; i++) r(i) = ref.e[family].vrr[i];

        RealVector final_direction;
        int edge;

        std::vector<RealVector> temp;

        curve_engine(ref.point, r, final_direction, temp, edge);
        if (temp.size() > 0){
            curves.push_back(temp);
            edges.push_back(edge);
        }

        temp.clear();
        curve_engine(ref.point, -r, final_direction, temp, edge);
        if (temp.size() > 0){
            curves.push_back(temp);
            edges.push_back(edge);
        }
    }

    return HUGONIOTCONTINUATION_CURVE_OK; 
}

//int HugoniotContinuation::non_local_branches(std::vector<Curve> &vc){
//    // Search method.
//    //
//    Hugoniot_Curve hc(flux, accum);

//    bool found = true;

//    RealVector initial_point = ref.point;

//    while (found){
//        std::vector< std::vector<RealVector> > local_hugoniot;
//        curve(initial_point, local_hugoniot);

//        vc.push_back(local_hugoniot);

//        disable_grid(local_hugoniot, grid);

//        // Find the first segment using the search method.
//        //
//        std::vector<RealVector> search_hugoniot_curve;
//        hc.find_segment(ref, grid, search_hugoniot_curve);

//        if (search_hugoniot_curve.size() == 0) found = false;
//        else {
//            initial_point = search_hugoniot_curve[0];
//        }
//    }

//    return HUGONIOTCONTINUATION_CURVE_OK;
//}

// TODO: This method will migrate to Stone-something.
//
//void HugoniotContinuation::hugoniot_on_horizontal_side(double p, JetMatrix &j){
//    RealVector pvert(2);
//    pvert(0) = p;
//    pvert(1) = 0.0;

//    RealVector direction(2);
//    direction(0) = 1.0;
//    direction(1) = 0.0;

//    // Common stuff.
//}

//void HugoniotContinuation::hugoniot_on_hypotenuse(double p, JetMatrix &j){
//    RealVector pvert(2);
//    pvert(0) = p;
//    pvert(1) = 1.0 - p;

//    RealVector direction(2);
//    direction(0) =  1.0;
//    direction(1) = -1.0;

//    // Common stuff.
//}

//
//
//

#include "Three_Phase_Boundary.h"
void HugoniotContinuation::hugoniot_on_side(double p, int side, double &RHv, double &RHv_prime){
    RealVector point(2);
    RealVector direction(2);

    if (side == THREE_PHASE_BOUNDARY_SW_ZERO){
        // Vertical.
        //
        point(0) = 0.0;
        point(1) = p;

        direction(0) = 0.0;
        direction(1) = 1.0;
    }
    else if (side == THREE_PHASE_BOUNDARY_SO_ZERO){
        // Horizontal.
        //
        point(0) = p;
        point(1) = 0.0;

        direction(0) = 1.0;
        direction(1) = 0.0;
    }
    else if (side == THREE_PHASE_BOUNDARY_SG_ZERO){
        // Diagonal.
        //
        point(0) = p;
        point(1) = 1.0 - p;

        direction(0) =  1.0;
        direction(1) = -1.0;
    }

    JetMatrix F_jet(point.size()), G_jet(point.size());
    f->jet(point, F_jet, 1);
    g->jet(point, G_jet, 1);

    RealVector H;
    DoubleMatrix nablaH;
    jet_Hugoniot(F_jet.function(), F_jet.Jacobian(), 
                 G_jet.function(), G_jet.Jacobian(),
                 H, nablaH);

    RHv = H(0);

    // nablaH comes transposed out of jet_Hugoniot. 
    RHv_prime = (direction*nablaH)(0);

    return;
}

int HugoniotContinuation::Newton_on_side(double init_p, int side, double &p){
    p = init_p;
    double RHv, RHv_prime;

    int max_it = 100;
    int it = 0;

    double epsilon = 1e-10;
    double old_p;

    do {
        old_p = p;

        hugoniot_on_side(p, side, RHv, RHv_prime);
        p = p - RHv/RHv_prime;

        it++;
    } while (std::abs(old_p - p) > epsilon && it < max_it);

//    std::cout << "Newton_on_side. Iterations = " << it << std::endl;

    if (it < max_it) return HUGONIOTCONTINUATION_NEWTON_OK;
    else             return HUGONIOTCONTINUATION_NEWTON_ERROR;
}

bool HugoniotContinuation::find_a_point_on_a_side(int side, RealVector &p){
//    // See if this exists, and if the new Boundary deals with this as of now.
//    b->fill_side_points(side, side_points);

    RealVector p0(2), p1(2);
    int main_index;

    if (side == THREE_PHASE_BOUNDARY_SW_ZERO){
        p0(0) = 0.0;
        p0(1) = 0.0;

        p1(0) = 0.0;
        p1(1) = 1.0;

        main_index = 1;
    }
    else if (side == THREE_PHASE_BOUNDARY_SO_ZERO){
        p0(0) = 0.0;
        p0(1) = 0.0;

        p1(0) = 1.0;
        p1(1) = 0.0;

        main_index = 0;
    }
    else {
        p1(0) = 0.0;
        p1(1) = 1.0;

        p0(0) = 1.0;
        p0(1) = 0.0;

        main_index = 0;
    }

    int n = 20;
    double delta = 1.0/(double)(n - 1);

    std::vector<RealVector> side_points;
    for (int i = 0; i < n; i++){
        double alpha = ((double)i)*delta;
        side_points.push_back(alpha*p0 + (1.0 - alpha)*p1);
    }

    std::vector<double> side_p;
    for (int i = 0; i < n; i++) side_p.push_back(side_points[i](main_index));

    double Hugoniot_function;
    double der_Hugoniot_function; // Not used, pass a flag to tell that it does not need to be computed.

//    hugoniot_on_side(side_points[0], side, Hugoniot_function, der_Hugoniot_function);
    hugoniot_on_side(side_p[0], side, Hugoniot_function, der_Hugoniot_function);

    for (int i = 1; i < side_points.size(); i++){
        double next_Hugoniot_function;

//        hugoniot_on_side(side_points[i], side, next_Hugoniot_function, der_Hugoniot_function);
        hugoniot_on_side(side_p[i], side, next_Hugoniot_function, der_Hugoniot_function);

        if (Hugoniot_function*next_Hugoniot_function < 0.0){
            double init_p = .5*(side_p[i - 1] + side_p[i]);

            double p_seed;
            Newton_on_side(init_p, side, p_seed);

            p.resize(2);
            if (side == THREE_PHASE_BOUNDARY_SW_ZERO){
                p(0) = 0.0;
                p(1) = p_seed;
            }
            else if (side == THREE_PHASE_BOUNDARY_SO_ZERO){
                p(0) = p_seed;
                p(1) = 0.0;
            }
            else{
                p(0) = p_seed;
                p(1) = 1.0 - p_seed;
            }

            return true;
        }
        else {
            Hugoniot_function = next_Hugoniot_function;
        }
    }

    return false;
}

bool HugoniotContinuation::find_a_point_on_a_side(int side, 
                                                  const std::vector<RealVector> &side_points,
                                                  const std::vector<bool> &use_this_segment,
                                                  RealVector &p){
//    // See if this exists, and if the new Boundary deals with this as of now.
//    b->fill_side_points(side, side_points);

    RealVector p0(2), p1(2);
    int main_index;

    if (side == THREE_PHASE_BOUNDARY_SW_ZERO){
        p0(0) = 0.0;
        p0(1) = 0.0;

        p1(0) = 0.0;
        p1(1) = 1.0;

        main_index = 1;
    }
    else if (side == THREE_PHASE_BOUNDARY_SO_ZERO){
        p0(0) = 0.0;
        p0(1) = 0.0;

        p1(0) = 1.0;
        p1(1) = 0.0;

        main_index = 0;
    }
    else {
        p1(0) = 0.0;
        p1(1) = 1.0;

        p0(0) = 1.0;
        p0(1) = 0.0;

        main_index = 0;
    }

//    int n = 20;
//    double delta = 1.0/(double)(n - 1);

//    std::vector<RealVector> side_points;
//    for (int i = 0; i < n; i++){
//        double alpha = ((double)i)*delta;
//        side_points.push_back(alpha*p0 + (1.0 - alpha)*p1);
//    }

    int n = side_points.size();

    std::vector<double> side_p;
    for (int i = 0; i < n; i++) side_p.push_back(side_points[i](main_index));

    double Hugoniot_function;
    double der_Hugoniot_function; // Not used, pass a flag to tell that it does not need to be computed.

//    hugoniot_on_side(side_points[0], side, Hugoniot_function, der_Hugoniot_function);
    hugoniot_on_side(side_p[0], side, Hugoniot_function, der_Hugoniot_function);

    for (int i = 1; i < side_points.size(); i++){
        if (!use_this_segment[i - 1]) continue;

        double next_Hugoniot_function;

//        hugoniot_on_side(side_points[i], side, next_Hugoniot_function, der_Hugoniot_function);
        hugoniot_on_side(side_p[i], side, next_Hugoniot_function, der_Hugoniot_function);

        if (Hugoniot_function*next_Hugoniot_function < 0.0){
            double init_p = .5*(side_p[i - 1] + side_p[i]);

            double p_seed;
            Newton_on_side(init_p, side, p_seed);

            p.resize(2);
            if (side == THREE_PHASE_BOUNDARY_SW_ZERO){
                p(0) = 0.0;
                p(1) = p_seed;
            }
            else if (side == THREE_PHASE_BOUNDARY_SO_ZERO){
                p(0) = p_seed;
                p(1) = 0.0;
            }
            else{
                p(0) = p_seed;
                p(1) = 1.0 - p_seed;
            }

            return true;
        }
        else {
            Hugoniot_function = next_Hugoniot_function;
        }
    }

    return false;
}

//// The lines below belong to a block, should be commented/uncommented as a whole.

//bool HugoniotContinuation::point_on_side_segment(std::vector<RealVector> &point_on_boundary, int side, const RealVector &segment_begin, const RealVector &segment_end){

//    std::vector<RealVector> temp_point_on_boundary;
//    for (int i = 0; i < point_on_boundary.size(); i++){
//        if (!boundary->point_on_side_segment(point_on_boundary[i], side, segment_begin, segment_end)) temp_point_on_boundary.push_back(point_on_boundary[i]);
//    }

//    point_on_boundary = temp_point_on_boundary;

//    return ;
//}

//int HugoniotContinuation::points_on_side(std::vector<RealVector> &point_on_boundary, int side, int n, std::vector<RealVector> &seeds){
//    std::vector<RealVector> side_points;

//    boundary->fill_side_points(side, n, side_points);

//    for (int i = 0; i < side_points.size() - 1; i++){
//        
//    }

//    return;
//}

//void HugoniotContinuation::identify_segments(const std::vector<Curve> &existing_curves, 
//                                             std::vector<std::vector<bool> > check_this_segment){

//    check_this_segment.clear();
//    // Reserve space here for check_this_segment.


//    for (int i = 0; i < existing_curves.sides(); i++){
//        // Check, update check_this_segment.
//    }

//    return;
//}


//int HugoniotContinuation::(const std::vector<Curve> &existing_curves, 
//                           std::vector<std::vector<bool> > check_this_segment,
//                           int max_number_of_curves,
//                           std::vector<Curve> &curves, std::vector<RealVector> &seeds){


//    int number_of_curves = 0;

//    // Obtain the segments that form each of the sides of the boundary.
//    //
//    std::vector<std::vector<RealVector> > side_segments_vertices;
//    boundary->side_segments_vertices(side_segments_vertices); // TODO: Create this method.

//    int number_of_sides = boundary->side_segments_vertices.size();

//    for (int i = 0; i < number_of_sides; i++){
//        // Number of segments in each side.
//        int number_of_vertices = side_segments_vertices[i].size();
//        std::vector<RealVector> &this_side_vertices = side_vertices_segments[i];

//        // Fill the list with the value of the Rankine-Hugoniot expression in
//        // all the points that form this side.
//        //
//        std::vector<double> RH_exp(number_of_vertices);
//        #pragma omp parallel for schedule(dynamic)
//        for (int j = 0; j < number_of_vertices; j++){
//            RH_exp[j] = Rankine_Hugoniot_expression(this_side_vertices[j]);
//        }
//        
//        int number_of_segments = number_of_vertices - 1;
//        for (int j = 0; j < number_of_segments; j++){
//            if (check_this_segment[i][j]){
//                if (RH_exp[j]*RH_exp[j + 1] < 0.0) {
//                    check_this_segment[i][j] = false; // Avoid this segment in the future

//                    // Find the seed and the curve here.
//                    Curve non_local_curve;
//                    curve_engine(..., non_local_curve);

//                    curves.push_back(non_local_curve);

//                    // Add the beginning of new_curve as a seed.
//                    //
//                    seeds.push_back(non_local_curve.curve.front());

//                    // The cell where the last point of this curve lies must be marked as invalid.
//                    identify_segments(const std::vector<Curve> &existing_curves, 
//                                             std::vector<std::vector<bool> > check_this_segment);

//                    //
//                    number_of_curves++;
//                    if (number_of_curves >= max_number_of_curves) return;
//                }
//            } 
//            
//        }
//    }

//    return;
//}


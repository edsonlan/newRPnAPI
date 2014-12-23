#include "NewRarefactionCurve.h"

NewRarefactionCurve::NewRarefactionCurve(Equations *eq, const Boundary *bb){
    equations_ = eq;
    b = bb;
}

NewRarefactionCurve::NewRarefactionCurve(const AccumulationFunction *gg, const FluxFunction *ff, const Boundary *bb){
    f = ff;
    g = gg;
    b = bb;
}

NewRarefactionCurve::~NewRarefactionCurve(){
}

int NewRarefactionCurve::field(int *neq, double *xi, double *in, double *out, int *obj, double* /* Not used */){
    // Recover the object.
    //
    NewRarefactionCurve *rar = (NewRarefactionCurve*)obj;

    // Verify that the point given by the ODE solver is valid.
    //
    RealVector point(*neq, in);
    if (!rar->b->inside(point)){
        //std::cout << "Field. point = " << point << " is OUTSIDE!" << std::endl;

        return FIELD_POINT_OUTSIDE_DOMAIN;
    }

//    // Proceed. Find the eigenpairs.
//    //
//    JetMatrix F_jet(*neq), G_jet(*neq);
//    rar->f->jet(point, F_jet, 1);
//    rar->g->jet(point, G_jet, 1);


    DoubleMatrix DF, DG;
    rar->equations_->constrained_conservation_laws_jet(point, DF, DG);

    std::vector<eigenpair> e;
    Eigen::eig(*neq, DF.data(), DG.data(), e);
 
    // TODO: Verify that the eigenvalue is not complex.
  
    // Extract right-eigenvector, verify that it points in the correct direction.
    //
    RealVector rm(*neq, e[rar->family].vrr.data());
    if (rm*rar->reference_vector < 0.0) rm = -rm;

    // Output.
    //
    for (int i = 0; i < *neq; i++) out[i] = rm(i);

    return FIELD_OK;
}

int NewRarefactionCurve::Jacobian_field(){

//subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
//c   dimension y(neq), pd(nrowpd,neq)
//c which supplies df/dy by loading pd as follows..
//c for a full jacobian (mf = 21), load pd(i,j) with df(i)/dy(j),
//c the partial derivative of f(i) with respect to y(j).  (ignore the
//c ml and mu arguments in this case.)

//int NewRarefactionCurve::field_Jacobian(, int *object){
//    // Recover the object.
//    //
//    NewRarefactionCurve *rar = (NewRarefactionCurve*)obj;

//    // Verify that the point given by the ODE solver is valid.
//    //
//    RealVector point(*neq, in);
//    if (!rar->b->inside(point)){
//        std::cout << "NewRarefactionCurve\'s field\'s Jacobian. point = " << point << " is OUTSIDE!" << std::endl;

//        return JACOBIAN_FIELD_POINT_OUTSIDE_DOMAIN;
//    }

//    

//    // Output.
//    //
//    for (int i = 0; i < *neq; i++){
//        for (int j = 0; j < *nrowpd; j++){
//            pd[i*] = J(j, i);
//        }
//    }

//    return JACOBIAN_FIELD_OK;
//}

}

// TODO: Remember this method will ultimately be moved to another class.
//
// TODO: Is the family really necessary?
//
double NewRarefactionCurve::directional_derivative(const RealVector &p, int fam, const RealVector &ref){
    int n = p.size();

    JetMatrix F_jet(n), G_jet(n);
    f->jet(p, F_jet, 2);
    g->jet(p, G_jet, 2);    

    DoubleMatrix F_J = F_jet.Jacobian();
    DoubleMatrix G_J = G_jet.Jacobian();

    std::vector<DoubleMatrix> F_H = F_jet.Hessian();
    std::vector<DoubleMatrix> G_H = G_jet.Hessian();

    std::vector<eigenpair> e; 
    Eigen::eig(n, F_J.data(), G_J.data(), e);

    double lambda = e[fam].r; 
    RealVector rm(n, e[fam].vrr.data()); // Right eigenvector
    RealVector lm(n, e[fam].vlr.data()); // Left  eigenvector

    // Verify that rm points in the right direction.
    // Since rm will be used below, it is necessary that this change is permanent.
    //
    if (rm*ref < 0.0) rm = -rm;

    RealVector h(F_H.size());

    for (int i = 0; i < F_H.size(); i++) h(i) = rm*((F_H[i] - lambda*G_H[i])*rm);
    
    return (lm*h)/(lm*(G_J*rm));
}

int NewRarefactionCurve::initialize(const RealVector &p, int fam, const RealVector &direction, RealVector &ref, double &dd){
    ref = direction;

    dd = directional_derivative(p, fam, ref);

    return RAREFACTION_INIT_OK;
}

int NewRarefactionCurve::initialize(const RealVector &p, int fam, int increase, RealVector &ref, double &dd){
    // If the directional derivative turns out to be null, the test below should be used.
    // In other words, when the directional derivative is operational again, add this test.
    // This will happen when at an inflection. Two curves or none shall be computed from the point.
    //
    DoubleMatrix DF, DG;
    equations_->constrained_conservation_laws_jet(p, DF, DG);
    
    std::vector<eigenpair> e;
    Eigen::eig(p.size(), DF.data(), DG.data(), e);

    double epsilon = 1e-6;
    RealVector r(p.size(), e[fam].vrr.data());
    double lambda = e[fam].r;

    RealVector pp = p + epsilon*r;
    double lp;
    {
        DoubleMatrix DFp, DGp;
        equations_->constrained_conservation_laws_jet(pp, DFp, DGp);
    
        std::vector<eigenpair> e;
        Eigen::eig(pp.size(), DFp.data(), DGp.data(), e);

        lp = e[fam].r;
    }

    RealVector pm = p - epsilon*r;
    double lm;
    {
        DoubleMatrix DFm, DGm;
        equations_->constrained_conservation_laws_jet(pm, DFm, DGm);
    
        std::vector<eigenpair> e;
        Eigen::eig(pm.size(), DFm.data(), DGm.data(), e);

        lm = e[fam].r;
    }

    if (lp > lambda && lambda > lm){
        if (increase == RAREFACTION_SPEED_SHOULD_INCREASE) ref = r;
        else                                               ref = -r;
    }
    else if (lp < lambda && lambda < lm){
        if (increase == RAREFACTION_SPEED_SHOULD_INCREASE) ref = -r;
        else                                               ref = r;
    }
    else return RAREFACTION_INIT_ERROR;

    //dd = directional_derivative(p, fam, ref);
    dd = 0.0; // Until the directional derivative is active again.

    return RAREFACTION_INIT_OK;
}

//// TODO: Errors can occur when initializing.
////       Deal with it.
//int NewRarefactionCurve::initialize(const RealVector &p, int fam, int increase, RealVector &ref, double &dd){
//    int n = p.size();

//    JetMatrix F_jet(n), G_jet(n);
//    f->jet(p, F_jet, 2);
//    g->jet(p, G_jet, 2);    

//    DoubleMatrix F_J = F_jet.Jacobian();
//    DoubleMatrix G_J = G_jet.Jacobian();

//    std::vector<DoubleMatrix> F_H = F_jet.Hessian();
//    std::vector<DoubleMatrix> G_H = G_jet.Hessian();

//    std::vector<eigenpair> e;
//    Eigen::eig(n, F_J.data(), G_J.data(), e);

//    double lambda = e[fam].r;
//    RealVector rm(n, e[fam].vrr.data()); // Right eigenvector
//    RealVector lm(n, e[fam].vlr.data()); // Left  eigenvector

//    RealVector h(F_H.size());

//    for (int i = 0; i < F_H.size(); i++) h(i) = rm*((F_H[i] - lambda*G_H[i])*rm);
//    
//    dd = (lm*h)/(lm*(G_J*rm));

//    if (((increase == RAREFACTION_SPEED_SHOULD_INCREASE) && (dd > 0.0)) ||
//        ((increase == RAREFACTION_SPEED_SHOULD_DECREASE) && (dd < 0.0))){
//        ref = rm;
//    }
//    else {
//        ref = -rm;
//        dd  = -dd;
//    }

//    return RAREFACTION_INIT_OK;
//}

int NewRarefactionCurve::inflection_signal_event(const RealVector & where, double & directional_derivative_measure, int *signal_object, int * /* reference_direction */){
    NewRarefactionCurve *rar = (NewRarefactionCurve*)signal_object;

    directional_derivative_measure = rar->directional_derivative(where, rar->family, rar->reference_vector);

    return BISECTION_FUNCTION_OK;
}

int NewRarefactionCurve::elliptic_region_signal_event_2D2D(const RealVector & where, double &discriminant, int *signal_object, int * /* reference_direction */){
    NewRarefactionCurve *rar = (NewRarefactionCurve*)signal_object;

    JetMatrix F_jet(2), G_jet(2);
    rar->f->jet(where, F_jet, 1);
    rar->g->jet(where, G_jet, 1);

    DoubleMatrix JF = F_jet.Jacobian();
    DoubleMatrix JG = G_jet.Jacobian();

    double A = JG(0, 0)*JG(1, 1) - JG(0, 1)*JG(1, 0);
    double B = (JG(0, 1)*JF(1, 0) + JF(0, 1)*JG(1, 0)) - (JG(0, 0)*JF(1, 1) + JF(0, 0)*JG(1, 1));
    double C = JF(0, 0)*JF(1, 1) - JF(0, 1)*JF(1, 0);

    discriminant = B*B - 4.0*A*C;

    return BISECTION_FUNCTION_OK;
}

// 3D2D
int NewRarefactionCurve::elliptic_region_signal_event_3D2D(const RealVector & where, double &discriminant, int *signal_object, int * /* reference_direction */){
    return BISECTION_FUNCTION_ERROR;

    NewRarefactionCurve *rar = (NewRarefactionCurve*)signal_object;

    JetMatrix F_jet(2), G_jet(2);
    rar->f->jet(where, F_jet, 1);
    rar->g->jet(where, G_jet, 1);

    DoubleMatrix JF = F_jet.Jacobian();
    DoubleMatrix JG = G_jet.Jacobian();

    double A = JG(0, 0)*JF(2, 2)*JG(1, 1) - JG(0, 0)*JF(1, 2)*JG(2, 1) + JG(1, 0)*JF(0, 2)*JG(2, 1) - JG(1, 0)*JF(2, 2)*JG(0, 1) + JG(2, 0)*JF(1, 2)*JG(0, 1) - JG(2, 0)*JF(0, 2)*JG(1, 1);

    double B = JF(0, 0)*JF(2, 2)*JG(1, 1) + JF(1, 2)*JF(0, 1)*JG(2, 0) + JF(2, 2)*JF(1, 1)*JG(0, 0) + JF(1, 0)*JF(0, 2)*JG(2, 1) - JF(1, 2)*JF(2, 1)*JG(0, 0) - JF(1, 0)*JF(2, 2)*JG(0, 1) +
               JF(0, 2)*JF(2, 1)*JG(1, 0) + JF(2, 0)*JF(1, 2)*JG(0, 1) - JF(2, 2)*JF(0, 1)*JG(1, 0) - JF(2, 0)*JF(0, 2)*JG(1, 1) - JF(0, 0)*JF(1, 2)*JG(2, 1) - JF(0, 2)*JF(1, 1)*JG(2, 0);

    double C = JF(0, 0)*JF(2, 2)*JF(1, 1) - JF(2, 0)*JF(0, 2)*JF(1, 1) - JF(0, 0)*JF(1, 2)*JF(2, 1) + JF(2, 0)*JF(1, 2)*JF(0, 1) + JF(1, 0)*JF(0, 2)*JF(2, 1) - JF(1, 0)*JF(2, 2)*JF(0, 1);

    discriminant = B*B - 4.0*A*C;

    return BISECTION_FUNCTION_OK;
}

void NewRarefactionCurve::all_eigenvalues(const RealVector &p, int fam, std::vector<std::complex<double> > &lambda){
//    std::vector<eigenpair> e;

//    Eigen::eig(p, f, g, e);

    DoubleMatrix DF, DG;
    equations_->constrained_conservation_laws_jet(p, DF, DG);
    
    std::vector<eigenpair> e;
    Eigen::eig(p.size(), DF.data(), DG.data(), e);

//    lambda.resize(e.size());
    for (int i = 0; i < e.size(); i++) lambda.push_back(std::complex<double>(e[i].r, e[i].i));

    return;
}

void NewRarefactionCurve::all_eigenvalues(const RealVector &p, int fam, RealVector &lambda){
    std::vector<std::complex<double> > complex_lambda;

    all_eigenvalues(p, fam, complex_lambda);

    lambda.resize(complex_lambda.size());
    for (int i = 0; i < lambda.size(); i++) lambda(i) = complex_lambda[i].real();

    return;
}

void NewRarefactionCurve::add_point_to_curve(const RealVector &p, Curve &curve){
    curve.back_pointer.push_back(curve.back_pointer.size() - 1);

    curve.curve.push_back(p);

    RealVector point_eigenvalues;
    all_eigenvalues(p, family, point_eigenvalues);

    curve.eigenvalues.push_back(point_eigenvalues);

    curve.speed.push_back(point_eigenvalues(family));

    return;
}

// TODO: Implement RAREFACTION_AS_ENGINE_FOR_INTEGRAL_CURVE & fill the inflection points.
//
int NewRarefactionCurve::curve(const RealVector &initial_point,
                            int curve_family,
                            int increase,
                            int type_of_rarefaction, // For itself or as engine for integral curve.
                            int should_initialize,
                            const RealVector *direction,
                            const ODE_Solver *odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
                            double deltaxi,
                            Curve &rarcurve,
                            std::vector<RealVector> &inflection_points, // Will these survive/be added to the Curve class?
                            RealVector &final_direction,
                            int &reason_why, // Similar to Composite.
                            int &edge){

    // Verify that initial_point is inside the boundary
    if (!b->inside(initial_point)){
        std::cout << "NewRarefactionCurve: The initial point " << initial_point << " is outside the domain! Aborting..." << std::endl; 

        return RAREFACTION_ERROR;
    }

    // Clear the curve, add some info.
    //
    rarcurve.clear();
    rarcurve.type = RAREFACTION_CURVE;
    rarcurve.family = curve_family;

    inflection_points.clear();

    // Set some parameters.
    //
    family = curve_family;

    // Directional derivative.
    //
    double dirdrv, previous_dirdrv;

    // Initialize: find the initial direction and directional derivative.
    //
    int info_initialize;

    if (should_initialize == RAREFACTION_INITIALIZE){
        info_initialize = initialize(initial_point, family, increase, reference_vector, dirdrv);
    }
    else {
        info_initialize = initialize(initial_point, family, *direction, reference_vector, dirdrv);
    }

    if (info_initialize == RAREFACTION_INIT_ERROR){
        std::cout << "NewRarefactionCurve: The initialization failed. Aborting..." << std::endl; 

        return RAREFACTION_ERROR;
    }

    // Initialize previous_dirdrv.
    //
    previous_dirdrv = dirdrv;

    // Add the initial point.
    //
    add_point_to_curve(initial_point, rarcurve);

    // Independent parameter.
    //
    double xi = 0.0;
    double next_xi = xi + deltaxi;

    // The points.
    //
    RealVector point(initial_point);
    RealVector next_point;

    while (true){
        if (rarcurve.curve.size() > 8000) return RAREFACTION_OK;
        int info_odesolver = odesolver->integrate_step(&NewRarefactionCurve::field, 
                                                      (int*)this, 0 /*double *function_data*/,
                                                      xi,      point,
                                                      next_xi, next_point);

        if (info_odesolver == ODE_SOLVER_ERROR){
            std::cout << "NewRarefactionCurve: The solver failed to find the next point (Error = " << info_odesolver << "). Aborting..." << std::endl; 

            return RAREFACTION_ERROR;
        }

        if (info_odesolver == ODE_SOLVER_POINT_OUTSIDE_DOMAIN){
            std::cout << "NewRarefactionCurve: The solver passed a point outside the domain to the field (Error = " << info_odesolver << "). Aborting..." << std::endl; 

            return RAREFACTION_ERROR;
        }

        // Has the rarefacion curve reached the boundary?
        //
        RealVector intersection_point;
        int info_intersect = b->intersection(next_point, point, intersection_point, edge);

        if (info_intersect == BOUNDARY_INTERSECTION_FOUND){
            // Add the point, etc.
            //
            add_point_to_curve(intersection_point, rarcurve);

            // The final direction is not needed, but will be added only for completeness.
            //
            final_direction = intersection_point - next_point;
            normalize(final_direction);
            rarcurve.final_direction = final_direction;

            rarcurve.last_point = intersection_point;

            reason_why = RAREFACTION_REACHED_BOUNDARY;
            rarcurve.reason_to_stop = RAREFACTION_REACHED_BOUNDARY;

            return RAREFACTION_OK;
        }

//        // Has the rarefaction curve reached an inflection? 
//        // Update the directional derivative first.
//        //
//        dirdrv = directional_derivative(next_point, family, reference_vector);

//        #ifdef TEST
//        {
//            std::vector<RealVector> v;
//            v.push_back(next_point);

//            std::vector<std::string> s;
//            std::stringstream ss;
////            ss << dirdrv;
//            ss << rarcurve.curve.size() - 1;
//            s.push_back(ss.str());

////            Curve2D *rar_point_for_canvas = new Curve2D(v, 1.0, 0.0, 0.0, s, CURVE2D_MARKERS | CURVE2D_INDICES);

////            canvas->add(rar_point_for_canvas);
////            scroll->add(ss.str().c_str(), canvas, rar_point_for_canvas);

////            std::cout << "Index = " << rarcurve.curve.size() - 1 << std::endl;

//            std::vector<eigenpair> e;
//            Eigen::eig(next_point, f, g, e);

//            int n = next_point.size();

////            std::cout << "    dirdrv   = " << dirdrv << std::endl;
////            std::cout << "    lambda 0 = " << e[0].r << std::endl;
////            std::cout << "    lambda 1 = " << e[1].r << std::endl;
////            std::cout << "    eigenvector 0 = " << RealVector(n, e[0].vrr.data()) << std::endl;
////            std::cout << "    eigenvector 1 = " << RealVector(n, e[1].vrr.data()) << std::endl;
//        }
//        #endif

//        if (dirdrv*previous_dirdrv < 0.0){
//            double bisection_epsilon = 1e-10; // Must be relative to the flux.

//            // Output here:
//            double c_t;
//            RealVector p_c;

//            // Invoke Bisection
//            int info_bisection = Bisection::bisection_method(xi,      point,
//                                                             next_xi, next_point,
//                                                             bisection_epsilon, 
//                                                             c_t,     p_c,
//                                                             &field, (int*)this, (double*)0,
//                                                             odesolver, // Should a different solver be used here? Another object of the same class?
//                                                             &inflection_signal_event, (int*)this, (int*)0);

//            if (info_bisection == BISECTION_FUNCTION_ERROR){
//                std::cout << "An error was reported by the signal function when called by Bisection. Leaving..." << std::endl;

//                rarcurve.reason_to_stop = RAREFACTION_ERROR;
//                return RAREFACTION_ERROR;
//            }
//            else if (info_bisection == BISECTION_EQUAL_SIGN){
//                std::cout << "Bisection detected that the signal event function has the same sign in both points. Leaving..." << std::endl;

//                rarcurve.reason_to_stop = RAREFACTION_ERROR;
//                return RAREFACTION_ERROR;
//            }
//            else if (info_bisection == BISECTION_CONVERGENCE_ERROR){
//                std::cout << "Bisection did not converge when computing the inflection. Leaving..." << std::endl;

//                rarcurve.reason_to_stop = RAREFACTION_ERROR;
//                return RAREFACTION_ERROR;
//            }
//            else {
//                std::cout << "Bisection converged when computing the inflection." << std::endl;

//                add_point_to_curve(p_c, rarcurve);

//                if (type_of_rarefaction == RAREFACTION){
//                    reason_why = RAREFACTION_REACHED_INFLECTION;
//                    rarcurve.reason_to_stop = RAREFACTION_REACHED_INFLECTION;

//                    rarcurve.last_point = p_c;
//                    rarcurve.final_direction = next_point - point; //rarcurve.curve.back() - rarcurve.curve[rarcurve.curve.size() - 2];
//                    normalize(rarcurve.final_direction);

//                    final_direction = rarcurve.final_direction;

//                    return RAREFACTION_OK;
//                }
//                else if (type_of_rarefaction == INTEGRAL_CURVE){
//                    inflection_points.push_back(p_c);
//                    next_point = p_c;
//                    dirdrv = 0.0;
//                }
//            }
//        
//        }

//        // Has the rarefaction curve reached the coincidence? 
//        //           
//        RealVector point_eigenvalues;
//        all_eigenvalues(next_point, family, point_eigenvalues);

////        for (int i = 0; i < point_eigenvalues.size() - 1; i++){
//            if (std::abs(point_eigenvalues(family) - point_eigenvalues(family + 1)) < 1e-10){
//                add_point_to_curve(next_point, rarcurve);

//                rarcurve.last_point = next_point;
//                rarcurve.final_direction = next_point - rarcurve.curve.back();
//                normalize(rarcurve.final_direction);

//                final_direction = rarcurve.final_direction;

//                reason_why = RAREFACTION_REACHED_COINCIDENCE_CURVE;
//                rarcurve.reason_to_stop = RAREFACTION_REACHED_COINCIDENCE_CURVE;

//                return RAREFACTION_OK;
//            }
////        }

        // Prepare for the next iteration...
        //
        xi = next_xi;
        next_xi += deltaxi;

        previous_dirdrv = dirdrv;

        // TODO: Use right-eigenvector here (if the points are too close).
        reference_vector = next_point - point;

        point = next_point;        

        // ...and store the results.
        //
        add_point_to_curve(point, rarcurve);
    }
}

//int NewRarefactionCurve::curve_from_boundary(const RealVector &initial_point, int side, 
//                  int curve_family,
//                  int increase,
//                  int type_of_rarefaction, // For itself or as engine for integral curve.
//                  const ODE_Solver *odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
//                  double deltaxi,
//                  Curve &rarcurve,
//                  std::vector<RealVector> &inflection_points, // Will these survive/be added to the Curve class?
//                  RealVector &final_direction,
//                  int &reason_why, // Similar to Composite.
//                  int &edge){

//    // Points to the interior of the domain from side s.
//    //
//    RealVector to_interior = b->side_transverse_interior(initial_point, side);

//    std::cout << "to_interior = " << to_interior << std::endl;

//    // Initialize.
//    //
//    RealVector initial_direction;
//    double dd;

//    int info_initialize = initialize(initial_point, curve_family, increase, initial_direction, dd);

//    std::cout << "info_initialize = " << info_initialize << std::endl;
//    std::cout << "dd = " << dd << ", initial_direction = " << initial_direction << std::endl;

//    if (info_initialize == RAREFACTION_INIT_ERROR) return RAREFACTION_ERROR;

//    std::cout << "initial_direction*to_interior = " << initial_direction*to_interior << std::endl;

//    // The rarefaction will be computed only if it can be computed from the boundary towards the interior
//    // of the domain (according to the requested value of increase).
//    // 
//    if (initial_direction*to_interior < 0.0) return RAREFACTION_ERROR;

//    int info = curve(initial_point, curve_family, increase, 
//                     type_of_rarefaction, RAREFACTION_DONT_INITIALIZE, 
//                     &initial_direction, odesolver, deltaxi,
//                     rarcurve, inflection_points, final_direction, 
//                     reason_why, edge);

//    return info;
//}

int NewRarefactionCurve::curve_from_boundary(const RealVector &initial_point, int side, 
                  int curve_family,
                  int increase,
                  int type_of_rarefaction, // For itself or as engine for integral curve.
                  const ODE_Solver *odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
                  double deltaxi,
                  Curve &rarcurve,
                  std::vector<RealVector> &inflection_points, // Will these survive/be added to the Curve class?
                  RealVector &final_direction,
                  int &reason_why, // Similar to Composite.
                  int &edge){

    // Points to the interior of the domain from side s.
    //
    RealVector to_interior = b->side_transverse_interior(initial_point, side);

    std::cout << "to_interior = " << to_interior << std::endl;

    // Find a point inside the domain, close to the initial point.
    //
    RealVector inner_point = initial_point + deltaxi*to_interior;

    // Find the lambdas.
    //
    std::vector<eigenpair> inner_e, initial_e;

    Eigen::fill_eigenpairs(f, g, inner_point, inner_e);
    Eigen::fill_eigenpairs(f, g, initial_point, initial_e);

    // Compare the lambdas.
    //
    double inner_lambda   = inner_e[curve_family].r;
    double initial_lambda = initial_e[curve_family].r;

    int n = initial_point.size();

    // The rarefaction will be computed only if it can be computed from the boundary towards the interior
    // of the domain (according to the requested value of increase).
    // 
    RealVector initial_direction;

    if ((inner_lambda > initial_lambda && increase == RAREFACTION_SPEED_SHOULD_INCREASE) ||
        (inner_lambda < initial_lambda && increase == RAREFACTION_SPEED_SHOULD_DECREASE)){

        RealVector rm = RealVector(n, initial_e[curve_family].vrr.data());

        initial_direction = (rm*to_interior > 0.0) ? rm : -rm; 

        RealVector temp;
        double dd;

        initialize(inner_point, curve_family, increase, temp, dd);
        std::cout << "dd = " << dd << ", temp = " << temp << std::endl;

        std::cout << "initial_direction = " << initial_direction << std::endl;
    }
    else return RAREFACTION_ERROR;

    int info = curve(inner_point, curve_family, increase, 
                     type_of_rarefaction, RAREFACTION_DONT_INITIALIZE, 
                     &initial_direction, odesolver, deltaxi,
                     rarcurve, inflection_points, final_direction, 
                     reason_why, edge);

    return info;
}


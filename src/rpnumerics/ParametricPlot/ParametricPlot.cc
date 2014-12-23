#include "ParametricPlot.h"


//
int ParametricPlot::find_initial_point_within_domain(RealVector (*f)(void*, double), void *obj, double &phi, double phi_final, double delta_phi, const Boundary *b, Curve &curve){
    RealVector point = (*f)(obj, phi);

    if (b->inside(point)){
        curve.curve.push_back(point);
        return INITIAL_POINT_FOUND;
    }
    else {
        RealVector old_point;

        while (phi <= phi_final){
           
            old_point = point;

            point = (*f)(obj, phi + delta_phi);

            RealVector boundary_point;
            int edge;

            int info = b->intersection(point, old_point, boundary_point, edge);

            if (info == BOUNDARY_INTERSECTION_FOUND){
                curve.curve.push_back(boundary_point);

//                std::cout << "Point at boundary: " << boundary_point << std::endl;
//                TestTools::pause("Point at boundary found. Check console.");

                double temp_phi;
                intersection(f, obj, b, 1.0 /*max_distance*/,
                             phi, phi + delta_phi, 
                             temp_phi);

                phi = temp_phi;

//                std::cout << "phi @ boundary = " << phi*180.0/M_PI << std::endl;
//                curve.curve.push_back(point);

                return INITIAL_POINT_FOUND;
            }

            phi += delta_phi;
        }
    }

    return INITIAL_POINT_NOT_FOUND;
}

void ParametricPlot::find_curve(RealVector (*f)(void*, double), 
                                bool (*f_asymptote)(void*, const RealVector&, const RealVector&), 
                                void *obj, 
                                double &phi, double phi_final, double max_distance, double delta_phi, const Boundary *b, Curve &curve){

    if (curve.curve.size() == 0) return;

    RealVector old_point = curve.curve.back();
//    std::cout << "Initial old_point = " << old_point << std::endl;

    int max_halvings = 20;
    int halvings = 0;

    // To increase.
    //
    int max_steps_current_delta_phi = 5;
    int steps_current_delta_phi = 0;

    while (phi <= phi_final){
        RealVector point = (*f)(obj, phi + delta_phi);

        if (norm(point - old_point) > max_distance*.01){
            while (norm(point - old_point) > max_distance*.01 && halvings < max_halvings){
                delta_phi *= .5;
                halvings++;

                point = (*f)(obj, phi + delta_phi);
            }
        }
        else {
            if (steps_current_delta_phi == max_steps_current_delta_phi){
                if (curve.curve.size() > 2){
                    int last_position = curve.curve.size() - 1;

                    RealVector old_direction = curve.curve[last_position - 1] - curve.curve[last_position - 2];
                    normalize(old_direction);

                    RealVector     direction = curve.curve[last_position] - curve.curve[last_position - 1];
                    normalize(direction);

                    if (old_direction*direction > MAX_COS_ANGLE){
                        steps_current_delta_phi = 0;
                        delta_phi *= SQRT_TWO;
                    }
                }
            }
            else {
                steps_current_delta_phi++;
            }

            point = (*f)(obj, phi + delta_phi);
        }

        // Update phi.
        //
        phi += delta_phi;

        if (halvings >= max_halvings){
            if (f_asymptote != 0){
                if ((*f_asymptote)(obj, old_point, point)) return;
            }
        }

        if (b->inside(point)){
            curve.curve.push_back(point);

           
        }
        else {
            RealVector boundary_point;
            int edge;

            int info = b->intersection(point, old_point, boundary_point, edge);

//            std::cout << "Info = " << info << ", Segment " << point << "-" << old_point << " should intersect the boundary at " << boundary_point << std::endl;
//            std::cout << "    point is inside: " << b->inside(point) << ", old_point is inside = " << b->inside(old_point) << std::endl;
//            TestTools::pause("Intersection found.");

            // ParametricPlot sometimes interacts wrongly with RectBoundary, which sometimes returns a void intersection point.
            // I hope the if below solves this problem for all cases (it does for Quad2ExplicitHugoniotCurve). Another solution has to be found.
            // 
            // Background: Quad2ExplicitHugoniotCurve.
            //
            // Morante.
            if (info == BOUNDARY_INTERSECTION_FOUND) curve.curve.push_back(boundary_point);

            return;
        }

        old_point = point;
    }

    return;
}

void ParametricPlot::plot(RealVector (*f)(void*, double), bool (*f_asymptote)(void*, const RealVector&, const RealVector&),
                          void *obj, double phi_init, double phi_final, int n, const Boundary *b, std::vector<Curve> &curve){
    curve.clear();

    // Maximum possible distance within this Boundary.
    //
    double max_distance = b->max_distance();

    // Find the first point inside the boundary.
    //
    double delta_phi = (phi_final - phi_init)/(double)(n - 1);

    double phi = phi_init;

    while (phi <= phi_final){
        Curve temp;

        int info_find_initial_point = find_initial_point_within_domain(f, obj, phi, phi_final, delta_phi, b, temp);

        if (info_find_initial_point == INITIAL_POINT_FOUND) {
//            std::cout << "Initial point: " << temp.curve.back() << std::endl;
//            TestTools::pause("Initial point found, check console.");

            find_curve(f, f_asymptote, obj, phi, phi_final, max_distance, delta_phi, b, temp);

            if (temp.curve.size() > 0) curve.push_back(temp); // Should it be > 1?
        }
    }

    return;
}

void ParametricPlot::intersection(RealVector (*f)(void*, double), void *obj, const Boundary *b, double max_distance,
                                  double init_phi_p, double init_phi_q, 
                                  double &phi){
//    std::cout << "init_phi_p = " << init_phi_p << ", init_phi_q = " << init_phi_q << std::endl;
//    TestTools::pause();


//    if      ( inside(p) &&  inside(q)) return BOUNDARY_INTERSECTION_BOTH_INSIDE;
//    else if (!inside(p) && !inside(q)) return BOUNDARY_INTERSECTION_BOTH_OUTSIDE;
//    else {
//        int n = p.size();

        // Initialize the temporary parameters and their corresponding points.
        //
        double phi_p = init_phi_p;
        double phi_q = init_phi_q;

        RealVector pp = (*f)(obj, phi_p);
        RealVector qq = (*f)(obj, phi_q);

        // Switch the temporary points if need be, such that pp is inside and qq is outside.
        //
        if (!b->inside(pp)){
            RealVector temp(pp);
            pp = qq;
            qq = temp;

            double temp_phi = phi_p;
            phi_p = phi_q;
            phi_q = temp_phi;

//            std::cout << "Inverted." << std::endl;
        }

//        TestTools::pause("Cycle");

        // Iterate while the distance between the points is large.
        //
        int count = 0;

        while (norm(pp - qq) > max_distance*.001 && count < 100) {
            count++;

            phi = .5*(phi_p + phi_q);
            RealVector r = (*f)(obj, phi);

            if (b->inside(r)){
                pp = r;
                phi_p = phi;
            }
            else {
                qq = r;
                phi_q = phi;
            }
        }

        return;

//    }
}

//void ParametricPlot::segmented_curve(RealVector (*f)(void*, double), void *obj, const Boundary *b,
//                                     double phi_init, double phi_final,
//                                     std::vector<RealVector> &sc){
//    
//    

//    return;
//}


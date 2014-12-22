#include "Utilities.h"

void Utilities::points_to_segments(std::vector<RealVector> &v){
    std::vector<RealVector> temp;

    if (v.size() > 0){
        temp.push_back(v.front());

        for (int i = 1; i < v.size() - 1; i++){
            temp.push_back(v[i]);
            temp.push_back(v[i]);
        }

        temp.push_back(v.back());
    }

    v.clear();
    for (int i = 0; i < temp.size(); i++) v.push_back(temp[i]);

    return;
}

int Utilities::bisection_on_segment(void* obj, double (*function_for_bisection)(void*, const RealVector&), const RealVector &p, const RealVector &q, RealVector &r){
    double fp = (*function_for_bisection)(obj, p);
    double fq = (*function_for_bisection)(obj, q);

    if (fp*fq > 0.0) return BISECTIONONSEGMENT_ERROR;

    RealVector q_minus_p(q - p);

    double alphap = 0.0;
    double alphaq = 1.0;

    int count = 0, max_count = 50;

    double alpha;

    while (std::abs(alphap - alphaq) > 1e-4 && count < max_count){
        alpha = .5*(alphap + alphaq);

        double f = (*function_for_bisection)(obj, p + alpha*q_minus_p);

        if (f*fp > 0.0){
            alphap = alpha;
            fp = f;
        }
        else {
            alphaq = alpha;
            fq = f;
        }

        count++;
    }

    if (count < max_count){
        r = p + alpha*q_minus_p;
        return BISECTIONONSEGMENT_OK;
    }
    else {
        return BISECTIONONSEGMENT_ERROR;
    }
}

//void Utilities::pick_point_from_continuous_curve(const std::vector<RealVector> &curve, const RealVector &p, RealVector &closest_point){
//    if (curve.size() < 2) return;

//    // Distance to a segment and the position of the first point of the segment in curve.
//    //
//    int pos = 0;
//    double min_distance = distance_point_line_2D(p, curve[0], curve[1]);

//    for (int i = 1; i < curve.size() - 1; i++){
//        double distance = distance_point_line_2D(p, curve[i], curve[i + 1]);

//        if (distance < min_distance){
//            pos = i;
//            min_distance = distance;
//        }
//    }

//    closest_point = project_point_onto_line_2D(p, curve[pos], curve[pos + 1]);

//    return;
//}

void Utilities::alpha_convex_combination_projection(const RealVector &p0, const RealVector &p1, const RealVector &p, double &alpha, RealVector &proj){
    RealVector seg_vec = p1 - p0;
    alpha = ((p - p0)*seg_vec)/norm(seg_vec);

    proj = p0 + alpha*seg_vec;

    return;
}

void Utilities::pick_point_from_continuous_curve(const std::vector<RealVector> &curve, const RealVector &p, RealVector &closest_point){
    int closest_segment_index;
    pick_point_from_continuous_curve(curve, p, closest_point, closest_segment_index);

//    if (curve.size() < 2) return;

//    // Distance to a point and the position of the first point of the segment in curve.
//    //
//    int pos = 0;
//    double min_distance = norm(p - curve[0]);

//    for (int i = 1; i < curve.size(); i++){
//        double distance = norm(p - curve[i]);

//        if (distance < min_distance){
//            pos = i;
//            min_distance = distance;
//        }
//    }

//    // http://paulbourke.net/geometry/pointlineplane/
//    //
//    if (pos == 0){
//        double alpha;
//        RealVector proj;
//        
//        alpha_convex_combination_projection(curve[pos], curve[pos + 1], p, alpha, proj);

//        if (validate_alpha(alpha)) closest_point = proj;
//        else                       closest_point = curve[pos];
//    }
//    else if (pos == curve.size() - 1){
//        double alpha;
//        RealVector proj;
//        
//        alpha_convex_combination_projection(curve[pos], curve[pos - 1], p, alpha, proj);

//        if (validate_alpha(alpha)) closest_point = proj;
//        else                       closest_point = curve[pos];
//    }
//    else {
//        // Check both segments.
//        //
//        double alpha_next, alpha_prev;
//        RealVector proj_next, proj_prev;
//        
//        alpha_convex_combination_projection(curve[pos], curve[pos + 1], p, alpha_next, proj_next);
//        alpha_convex_combination_projection(curve[pos], curve[pos - 1], p, alpha_prev, proj_prev);

//        bool alpha_next_valid = validate_alpha(alpha_next);
//        bool alpha_prev_valid = validate_alpha(alpha_prev);

//        if (alpha_next_valid && alpha_prev_valid){
//            double distance_next = norm2_squared(proj_next - p);
//            double distance_prev = norm2_squared(proj_prev - p);

//            if (distance_next < distance_prev) closest_point = proj_next;
//            else                               closest_point = proj_prev;
//        }
//        else {
//            if      (alpha_next_valid && !alpha_prev_valid) closest_point = proj_next;
//            else if (alpha_prev_valid && !alpha_next_valid) closest_point = proj_prev;
//            else                                            closest_point = curve[pos];
//        }

//    }

//    return;
}

void Utilities::pick_point_from_continuous_curve(const std::vector<RealVector> &curve, const RealVector &p, 
                                                 RealVector &closest_point, int &closest_segment_index){
    if (curve.size() < 2) return;

    // Distance to a point and the position of the first point of the segment in curve.
    //
    int pos = 0;
    double min_distance = norm(p - curve[0]);

    for (int i = 1; i < curve.size(); i++){
        double distance = norm(p - curve[i]);

        if (distance < min_distance){
            pos = i;
            min_distance = distance;
        }
    }

    // http://paulbourke.net/geometry/pointlineplane/
    //
    if (pos == 0){
        double alpha;
        RealVector proj;
        
        alpha_convex_combination_projection(curve[pos], curve[pos + 1], p, alpha, proj);

        if (validate_alpha(alpha)) closest_point = proj;
        else                       closest_point = curve[pos];

        closest_segment_index = pos;
    }
    else if (pos == curve.size() - 1){
        double alpha;
        RealVector proj;
        
        alpha_convex_combination_projection(curve[pos], curve[pos - 1], p, alpha, proj);

        if (validate_alpha(alpha)) closest_point = proj;
        else                       closest_point = curve[pos];

        closest_segment_index = pos - 1;
    }
    else {
        // Check both segments.
        //
        double alpha_next, alpha_prev;
        RealVector proj_next, proj_prev;
        
        alpha_convex_combination_projection(curve[pos], curve[pos + 1], p, alpha_next, proj_next);
        alpha_convex_combination_projection(curve[pos], curve[pos - 1], p, alpha_prev, proj_prev);

        bool alpha_next_valid = validate_alpha(alpha_next);
        bool alpha_prev_valid = validate_alpha(alpha_prev);

        if (alpha_next_valid && alpha_prev_valid){
            double distance_next = norm2_squared(proj_next - p);
            double distance_prev = norm2_squared(proj_prev - p);

            if (distance_next < distance_prev) closest_point = proj_next;
            else                               closest_point = proj_prev;
        }
        else {
            if (alpha_next_valid && !alpha_prev_valid){
                closest_point = proj_next;
                closest_segment_index = pos;
            }
            else if (alpha_prev_valid && !alpha_next_valid){
                closest_point = proj_prev;
                closest_segment_index = pos - 1;
            }
            else {
                closest_point = curve[pos];
                closest_segment_index = pos;
            }
        }
    }

    return;
}

void Utilities::pick_point_from_segmented_curve(const std::vector<RealVector> &curve, const RealVector &p, RealVector &closest_point){
    if (curve.size() < 2) return;

    double min_distance = std::numeric_limits<double>::infinity();
    RealVector p0, p1;

    double alpha;
    RealVector proj;

    for (int i = 0; i < curve.size()/2; i++){
        alpha_convex_combination_projection(curve[2*i], curve[2*i + 1], p, alpha, proj);

        double distance;

        if (alpha < 0.0)      distance = norm2_squared(p - curve[2*i]);
        else if (alpha > 1.0) distance = norm2_squared(p - curve[2*i + 1]);
        else                  distance = norm2_squared(p - proj);

        if (min_distance > distance){
            min_distance = distance;
            p0 = curve[2*i];
            p1 = curve[2*i + 1];
        }
    }

    alpha_convex_combination_projection(p0, p1, p, alpha, proj);

    if      (alpha < 0.0) closest_point = p0;
    else if (alpha > 1.0) closest_point = p1;
    else                  closest_point = proj;

    return;
}

void Utilities::pick_point_from_segmented_curve(const std::vector<RealVector> &curve, const RealVector &p, RealVector &closest_point, int &index_p0, int &index_p1){
    if (curve.size() < 2) return;

    double min_distance = std::numeric_limits<double>::infinity();

    double alpha;
    RealVector proj;

    for (int i = 0; i < curve.size()/2; i++){
        alpha_convex_combination_projection(curve[2*i], curve[2*i + 1], p, alpha, proj);

        double distance;

        if (alpha < 0.0)      distance = norm2_squared(p - curve[2*i]);
        else if (alpha > 1.0) distance = norm2_squared(p - curve[2*i + 1]);
        else                  distance = norm2_squared(p - proj);

        if (min_distance > distance){
            min_distance = distance;

            index_p0 = 2*i;
            index_p1 = 2*i + 1;
        }
    }

    RealVector p0 = curve[index_p0];
    RealVector p1 = curve[index_p1];

    alpha_convex_combination_projection(p0, p1, p, alpha, proj);

    if      (alpha < 0.0) closest_point = p0;
    else if (alpha > 1.0) closest_point = p1;
    else                  closest_point = proj;

    return;
}

int Utilities::find_point_on_level_curve(void *obj, double (*function_for_bisection)(void*, const RealVector &), const RealVector &p, const RealVector &q0, RealVector &point_on_level_curve){
    RealVector q = q0;
    RealVector vec = q0 - p;

    // Initialize.
    //
    double fp = (*function_for_bisection)(obj, p);
    double fq = (*function_for_bisection)(obj, q);

    //
    //
    int max_dup = 5;
    int dup = 0;

    while (fp*fq > 0.0 && dup < max_dup){
        std::cout << "while. dup = " << dup << std::endl;

        q = p + 2.0*vec; // Was: q = 2.0*(p + vec); ERROR!!!

        fq = (*function_for_bisection)(obj, q);
        vec = q - p;
        dup++;
    }

    // Bisection proper.
    //
    if (dup < max_dup){
        int info = bisection_on_segment(obj, function_for_bisection, p, q, point_on_level_curve);

        if (info == BISECTIONONSEGMENT_OK) return UTILITIES_BISECTION_OK;
    }
    
    return UTILITIES_BISECTION_ERROR;
}

void Utilities::pick_point_from_wavecurve(const WaveCurve &wavecurve, const RealVector &p, 
                                          int &curve_index, int &segment_index_in_curve, RealVector &closest_point, double &speed){
    std::vector<double> distance;
    std::vector<int> segment_index;
    std::vector<RealVector> point;

    for (int i = 0; i < wavecurve.wavecurve.size(); i++){
        RealVector closest_point;
        int closest_segment_index;

        pick_point_from_continuous_curve(wavecurve.wavecurve[i].curve, p, closest_point, closest_segment_index);

        segment_index.push_back(closest_segment_index);
        point.push_back(closest_point);
        distance.push_back(norm2_squared(closest_point - p));
    }

    double min_distance = std::numeric_limits<double>::infinity();
    int pos = 0;

    for (int i = 0; i < distance.size(); i++){
        if (min_distance > distance[i]){
            min_distance = distance[i];
            pos = i;
        }
    }

    closest_point = point[pos];
    curve_index = pos;
    segment_index_in_curve = segment_index[pos];

    double alpha;
    RealVector proj;

    // Point.
    //
    RealVector p0 = wavecurve.wavecurve[pos].curve[segment_index_in_curve];
    RealVector p1 = wavecurve.wavecurve[pos].curve[segment_index_in_curve + 1];

    alpha_convex_combination_projection(p0, p1, closest_point, alpha, proj);

    // Speed.
    //
    double speed0 = wavecurve.wavecurve[pos].speed[segment_index_in_curve];
    double speed1 = wavecurve.wavecurve[pos].speed[segment_index_in_curve + 1];

    speed = alpha*speed0 + (1.0 - alpha)*speed1;

    return;
}

// Created specifically for the Injection_to_Side case (Freddie's).
//
void Utilities::pick_last_point_from_wavecurve(const WaveCurve &wavecurve, const RealVector &p, 
                                               int &curve_index, int &segment_index_in_curve, RealVector &closest_point, double &speed){

    double min_distance = std::numeric_limits<double>::infinity();

    for (int i = 0; i < wavecurve.wavecurve.size(); i++){
        double distance = norm(p - wavecurve.wavecurve[i].curve.back());

        if (distance < min_distance){
            min_distance = distance;

            curve_index = i;

            // The last point in the curve. It will be decremented later on.
            //
            segment_index_in_curve = wavecurve.wavecurve[i].curve.size() - 1;
        }
    }

    // At least one point will be rejected now (the last one).
    //
    closest_point = wavecurve.wavecurve[curve_index].curve[segment_index_in_curve];
    while (norm(wavecurve.wavecurve[curve_index].curve.back() - closest_point) < 1e-4){
        segment_index_in_curve--;
        closest_point = wavecurve.wavecurve[curve_index].curve[segment_index_in_curve];
        speed = wavecurve.wavecurve[curve_index].speed[segment_index_in_curve];
    }

    return;
}

void Utilities::regularly_sampled_segment(const RealVector &p, const RealVector &q, int n, Curve &curve){
    curve.curve.clear();

    double delta = 1.0/(double)(n - 1);

    for (int i = 0; i < n; i++) curve.curve.push_back((p - q)*(delta*(double)(i)) + q);

    return;
}

int Utilities::Bhaskara(double b, double c, double &x1, double &x2){
    double disc = b*b - 4.0*c;
    if (disc < 0.0){
        std::cout << "Disc. = " << disc << std::endl;
        return BHASKARA_COMPLEX_ROOTS;
    }

    double sqrt_disc = sqrt(disc);

    // Rewrite the lines below if there is a 
    // catastrophic cancellation due to c being
    // too small.
    //
    x1 = 0.5*(-b + sqrt_disc);
    x2 = 0.5*(-b - sqrt_disc);

    if (sqrt_disc < 1e-8*std::abs(b)) return BHASKARA_DOUBLE_ROOTS;
    else                              return BHASKARA_TWO_DIFFERENT_ROOTS;
}

int Utilities::Bhaskara(double a, double b, double c, double &x1, double &x2){

    double disc = b*b - 4.0*a*c;
    if (disc < 0.0){
        std::cout << "Disc. = " << disc << std::endl;
        return BHASKARA_COMPLEX_ROOTS;
    }

    double sqrt_disc = sqrt(disc);

    // Rewrite the lines below if there is a 
    // catastrophic cancellation due to c being
    // too small.
    //
    int den = .5*(1.0/a);

    x1 = den*(-b + sqrt_disc);
    x2 = den*(-b - sqrt_disc);

    if (sqrt_disc < 1e-8*std::abs(b)) return BHASKARA_DOUBLE_ROOTS;
    else                              return BHASKARA_TWO_DIFFERENT_ROOTS;
}


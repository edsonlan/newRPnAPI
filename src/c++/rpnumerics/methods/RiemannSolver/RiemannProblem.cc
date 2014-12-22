#include "RiemannProblem.h"

RiemannProblem::RiemannProblem(){
}

RiemannProblem::~RiemannProblem(){
}

void RiemannProblem::half_profile(const WaveCurve &wavecurve, int initial_curve, int initial_point, int family,
                                  std::vector<RealVector> &phase_state, std::vector<double> &speed){
    int curve = initial_curve;
    int point = initial_point;

    std::cout << "curve = " << curve << std::endl;
    std::cout << "point = " << point << std::endl;

    while (true){
        if ((curve == -1 && point == 0) ||
            (curve == -1 && point == -1)
           ){
            // End this profile if the beginning of the first curve was reached.
            //
            break;
        }

        const Curve &current_curve = wavecurve.wavecurve[curve];
        RealVector current_point = current_curve.curve[point];
        double current_speed = current_curve.speed[point];

        phase_state.push_back(current_point);
        speed.push_back(current_speed);

        std::cout << "Riemann profile. Adding point " << point << " of curve " << curve << " (" << current_curve.type << ")" << std::endl;
        std::cout << "    speed = " << current_speed << std::endl;

        if (current_curve.type == SHOCK_CURVE){
            std::cout << "Riemann profile. Shock curve. speed = " << current_speed << std::endl;

            // If the current curve is a shock, go to the reference point and end this profile.
            //
            phase_state.push_back(wavecurve.reference_point.point);
//            speed.push_back(wavecurve.reference_point.e[family].r);
            speed.push_back(current_speed);

            break;
        }
        else if (current_curve.type == COMPOSITE_CURVE){
            std::cout << "Riemann composite"  << std::endl;
            std::cout << "point = " << point << ", curve = " << curve << std::endl;
        }

        int old_point = point;
        std::cout << "old_point = " << old_point << std::endl;

        point = current_curve.back_pointer[point];
        std::cout << "Next point: " << point << std::endl;

        curve = current_curve.back_curve_pointer[old_point];
        std::cout << "Next curve: " << curve << std::endl;    
        std::cout << "    current_curve.back_curve_pointer[old_point] = " << current_curve.back_curve_pointer[old_point] << std::endl;
    }

    return;
}

void RiemannProblem::profile(const WaveCurve &wavecurve1, int curve1, int point1, int family1,
                             const WaveCurve &wavecurve2, int curve2, int point2, int family2,
                             std::vector<RealVector> &phase_state, std::vector<double> &speed){

    // Check that the families and increases of both wavecurves are compatible.
    //
    if (wavecurve1.wavecurve[curve1].speed[point1] > wavecurve2.wavecurve[curve2].speed[point2]) return;

    phase_state.clear();
    speed.clear();

    // From L to M. Results must be reversed before the next step.
    //
    std::vector<RealVector> temp_phase_state;
    std::vector<double> temp_speed;

//        for (int i = 0; i < wavecurve1.wavecurve.size(); i++){
//            if (wavecurve1.wavecurve[i].type == COMPOSITE_CURVE){
//                std::cout << "Composite " << i << std::endl;
//                for (int j = 0; j < wavecurve1.wavecurve[i].curve.size(); j++){
//                    std::cout << "    j = " << j << ", back curve = " << wavecurve1.wavecurve[i].back_curve_pointer[j] << std::endl;
//                }
//            }
//        }

//    TestTools::pause("Check all the backcurves!");

    half_profile(wavecurve1, curve1, point1, family1, temp_phase_state, temp_speed);
    
    for (int i = temp_phase_state.size() - 1; i >= 0; i--){
        phase_state.push_back(temp_phase_state[i]);
        speed.push_back(temp_speed[i]);
    }

    // There is a small problem here, the profile retreats in speed. Try to find the bug. 
    //   
    RealVector p0(wavecurve1.wavecurve[curve1].curve[point1]);
    RealVector p1(wavecurve1.wavecurve[curve1].curve[point1 + 1]);

    RealVector q0(wavecurve2.wavecurve[curve2].curve[point2]);
    RealVector q1(wavecurve2.wavecurve[curve2].curve[point2 + 1]);

    RealVector r;
    double alpha, beta;

    bool info_seg = segment_segment_intersection(p0, p1, q0, q1, r, alpha, beta);
    if (!info_seg) return; // This REALLY should not happen.

    // TODO: There seems to be no universally good solution here.
    // If the curve is a rarefaction, interpolate. If it is a composite or a shock, draw a shock.

    // TODO:

    // Add the middle point with the speed of the left.
    //
    phase_state.push_back(r);
    //speed.push_back(alpha*wavecurve1.wavecurve[curve1].speed[point1] + (1.0 - alpha)*wavecurve1.wavecurve[curve1].speed[point1 + 1]);
    speed.push_back(wavecurve1.wavecurve[curve1].speed[point1]);

    // Add the middle point with the speed of the right.
    //
    phase_state.push_back(r);
//    speed.push_back(beta*wavecurve2.wavecurve[curve2].speed[point2] + (1.0 - beta)*wavecurve2.wavecurve[curve2].speed[point2 + 1]);
    speed.push_back(wavecurve2.wavecurve[curve2].speed[point2]);

    // From M to R.
    //
    half_profile(wavecurve2, curve2, point2, family2, phase_state, speed);

    return;
}

void RiemannProblem::all_increase_profile(const WaveCurve &wavecurve1, int curve1, int point1, int family1,
                                          const WaveCurve &wavecurve2, int curve2, int point2, int family2,
                                          std::vector<RealVector> &phase_state, std::vector<double> &speed){

    // Check that the families and increases of both wavecurves are compatible.
    //
    if (wavecurve1.wavecurve[curve1].speed[point1] > wavecurve2.wavecurve[curve2].speed[point2]) return;

    phase_state.clear();
    speed.clear();

    // From L to M. Results must be reversed before the next step.
    //
    std::vector<RealVector> temp_phase_state;
    std::vector<double> temp_speed;

    half_profile(wavecurve1, curve1, point1, family1, temp_phase_state, temp_speed);
    
    for (int i = temp_phase_state.size() - 1; i >= 0; i--){
        phase_state.push_back(temp_phase_state[i]);
        speed.push_back(temp_speed[i]);
    }

//    // There is a small problem here, the profile retreats in speed. Try to find the bug. 
//    //   
//    RealVector p0(wavecurve1.wavecurve[curve1].curve[point1]);
//    RealVector p1(wavecurve1.wavecurve[curve1].curve[point1 + 1]);

//    RealVector q0(wavecurve2.wavecurve[curve2].curve[point2]);
//    RealVector q1(wavecurve2.wavecurve[curve2].curve[point2 + 1]);

//    RealVector r;
//    double alpha, beta;

//    bool info_seg = segment_segment_intersection(p0, p1, q0, q1, r, alpha, beta);
//    if (!info_seg) return; // This REALLY should not happen.

    // TODO: There seems to be no universally good solution here.
    // If the curve is a rarefaction, interpolate. If it is a composite or a shock, draw a shock.

    // TODO:

    // Add the middle point with the speed of the left.
    //
//    phase_state.push_back(r);
    //speed.push_back(alpha*wavecurve1.wavecurve[curve1].speed[point1] + (1.0 - alpha)*wavecurve1.wavecurve[curve1].speed[point1 + 1]);

    phase_state.push_back(wavecurve1.wavecurve[curve1].curve[point1]);
    speed.push_back(wavecurve1.wavecurve[curve1].speed[point1]);

    // Add the middle point with the speed of the right.
    //
//    phase_state.push_back(r);
////    speed.push_back(beta*wavecurve2.wavecurve[curve2].speed[point2] + (1.0 - beta)*wavecurve2.wavecurve[curve2].speed[point2 + 1]);
//    speed.push_back(wavecurve2.wavecurve[curve2].speed[point2]);

    phase_state.push_back(wavecurve2.wavecurve[0].curve[0]);
    speed.push_back(wavecurve2.wavecurve[0].speed[0]);

    // From M to R.
    //
    temp_phase_state.clear();
    temp_speed.clear();

    half_profile(wavecurve2, curve2, point2, family2, temp_phase_state, temp_speed);
    std::cout << "From R to ref. temp_phase_state.size() = " << temp_phase_state.size() << std::endl;


    for (int i = temp_phase_state.size() - 1; i >= 0; i--){
        phase_state.push_back(temp_phase_state[i]);
        speed.push_back(temp_speed[i]);
    }

    return;
}



#include "WaveCurve.h"

void WaveCurve::init(const WaveCurve &orig){
    family = orig.family;

    increase = orig.increase;

    reference_point = orig.reference_point;

    for (int i = 0; i < orig.wavecurve.size(); i++) wavecurve.push_back(orig.wavecurve[i]);

    return;
}

WaveCurve::WaveCurve(){
}

WaveCurve::WaveCurve(const WaveCurve &orig){
    init(orig);
}

WaveCurve::WaveCurve(const WaveCurve *orig){
    init(*orig);
}

WaveCurve::~WaveCurve(){
}

WaveCurve WaveCurve::operator=(const WaveCurve &orig){
    if (this != &orig) init(orig);

    return *this;
}

void WaveCurve::add(const Curve &c){
    std::cout << "WaveCurve, will try to push back." << std::endl;

    wavecurve.push_back(c);

    std::cout << "WaveCurve, successfully pushed back." << std::endl;

    return;
}

void WaveCurve::insert_rarefaction_composite_pair(int composite_curve_index, int composite_insertion_point_index, const RealVector &rar_point, const RealVector &cmp_point){
    if (wavecurve[composite_curve_index].type == COMPOSITE_CURVE){
        // Shift the indices on this composite and its associated rarefaction.

        if (composite_insertion_point_index >= 0 && composite_insertion_point_index < wavecurve[composite_curve_index].curve.size() - 1){
            // Index of the related rarefaction curve (to simplify the access).
            //
            int rarefaction_index = wavecurve[composite_curve_index].back_curve_index;
            int rarefaction_point_index = wavecurve[composite_curve_index].back_pointer[composite_insertion_point_index]; // See pic IMG_20140112_153146316.jpg. TODO: Explain.
            
            // Insert the point into the rarefaction curve.
            //
            std::vector<RealVector>::iterator rar_it = wavecurve[rarefaction_index].curve.begin();
            wavecurve[rarefaction_index].curve.insert(rar_it + rarefaction_point_index + 1, rar_point);

            // Update the back pointers for the rarefaction.
            //
            for (int i = rarefaction_point_index; i < wavecurve[rarefaction_index].curve.size(); i++) wavecurve[rarefaction_index].back_pointer[i] = i - 1; // TODO: Check if it's ok.

            // Insert the point into the composite curve.
            //
            std::vector<RealVector>::iterator cmp_it = wavecurve[composite_curve_index].curve.begin();
            wavecurve[composite_curve_index].curve.insert(cmp_it + composite_insertion_point_index, cmp_point);

            // Update the back pointers for the composite.
            //
            for (int i = composite_insertion_point_index; i >= 0; i--) wavecurve[composite_curve_index].back_pointer[i] = wavecurve[composite_curve_index].back_pointer[i - 1];

            // Find all the composite curves that also depend on this rarefaction and propagate the shifting of indices.
            // The search must be backwards (approaching the rarefaction), because the points in the rarefaction
            // that are shifted due to an insertion affect the composite curves that are closer to the rarefaction.
            //
            int pos = composite_curve_index - 1;
            while (pos > rarefaction_index){
                if (wavecurve[pos].type ==  COMPOSITE_CURVE){
                }

                pos--;
            }
        }
    }

    return;
}

void WaveCurve::insert_rarefaction_composite_pair(int composite_curve_index, int composite_insertion_point_index, const RealVector &rarcmp_point){
    int n = rarcmp_point.size()/2;

    insert_rarefaction_composite_pair(composite_curve_index, composite_insertion_point_index, RealVector(0, n, rarcmp_point), RealVector(n, n, rarcmp_point));

    return;
}

void WaveCurve::validate_compatibility(WaveCurve &other_wavecurve, int curve_index, int point_index, int direction){
    double d;
    if (direction == SPEED_INCREASE){
        d =  1.0;
    }
    else {
        d = -1.0;
    }

    double speed = wavecurve[curve_index].speed[point_index];

    for (int i = 0; i < other_wavecurve.wavecurve.size(); i++){
        other_wavecurve.wavecurve[i].compatible.clear();

        for (int j = 0; j < other_wavecurve.wavecurve[i].speed.size(); j++){
            if ((speed - other_wavecurve.wavecurve[i].speed[j])*d < 0.0){
                other_wavecurve.wavecurve[i].compatible.push_back(true);
            }
            else {
                other_wavecurve.wavecurve[i].compatible.push_back(false);
            }
        }
    }

    return;
}

void WaveCurve::clear(){
    for (int i = 0; i < wavecurve.size(); i++) wavecurve[i].clear();

    return;
}

// No cleaning of vectors here.
//
void WaveCurve::half_speed_map(int begin, int end, std::vector<double> &arc_length, std::vector<double> &speed, std::vector<RealVector> &eigenvalues) const {
    // No cleaning here.
    //
    double distance = 0.0;
    RealVector prev_point = wavecurve[begin].curve.front();

    for (int i = begin; i < end; i++){
        for (int j = 0; j < wavecurve[i].curve.size(); j++){
            distance += norm(wavecurve[i].curve[j] - prev_point);
            arc_length.push_back(distance);

            speed.push_back(wavecurve[i].speed[j]);

            eigenvalues.push_back(wavecurve[i].eigenvalues[j]);

            prev_point = wavecurve[i].curve[j];
        }
    }

    return;
}

void WaveCurve::speed_map(std::vector<double> &arc_length, std::vector<double> &speed, std::vector<RealVector> &eigenvalues, RealVector &reference_eigenvalues) const {
    arc_length.clear();
    speed.clear();
    eigenvalues.clear();

    // First the second half-wavecurve will be processed. The arc length will be deemed negative, starting from the
    // initial point.
    //
    std::vector<double> temp_arc_length, temp_speed;
    std::vector<RealVector> temp_eigenvalues;

    half_speed_map(beginnig_of_second_half, wavecurve.size(), temp_arc_length, temp_speed, temp_eigenvalues);

    // Revert (using negative arclengths)
    for (int i = temp_speed.size() - 1; i >= 0; i--){
        arc_length.push_back(-temp_arc_length[i]);
        speed.push_back(temp_speed[i]);
        eigenvalues.push_back(temp_eigenvalues[i]);
    }

    // Now compute the first half-wavecurve. Notice that half_speed_map will not clean the vectors
    // so it can be used like this.
    //
    half_speed_map(0, beginnig_of_second_half, arc_length, speed, eigenvalues);

    reference_eigenvalues.resize(reference_point.e.size());
    for (int i = 0; i < reference_point.e.size(); i++) reference_eigenvalues(i) = reference_point.e[i].r;

    return;
}

void WaveCurve::half_speed_map(int type, int begin, int end, std::vector<Curve> &arclength_speed, std::vector<Curve> &arclength_eigenvalues) const {
    // No cleaning here.
    //
    double distance = 0.0;
    RealVector prev_point = wavecurve[begin].curve.front();
    RealVector point(2);

    RealVector eigen(wavecurve[begin].eigenvalues.front().size() + 1);

    double factor = (type == WAVECURVE_POSITIVE_ARCLENGTH) ? 1.0 : -1.0;

    for (int i = begin; i < end; i++){
        Curve temp_arclength_speed, temp_arclength_eigenvalues;
        temp_arclength_speed.type = wavecurve[i].type;
        temp_arclength_eigenvalues.type = wavecurve[i].type;

        for (int j = 0; j < wavecurve[i].curve.size(); j++){
            // For all.
            //
            distance += norm(wavecurve[i].curve[j] - prev_point);

            // For arclength_speed.
            //
            point(0) = factor*distance;
            point(1) = wavecurve[i].speed[j];
            temp_arclength_speed.curve.push_back(point);

            // For the eigenvalues.
            //
            eigen(0) = factor*distance;
            for (int k = 0; k < wavecurve[i].eigenvalues[j].size(); k++) eigen(1 + k) = wavecurve[i].eigenvalues[j](k);
            temp_arclength_eigenvalues.curve.push_back(eigen);
          
            // Update.
            //
            prev_point = wavecurve[i].curve[j];
        }

        arclength_speed.push_back(temp_arclength_speed);
        arclength_eigenvalues.push_back(temp_arclength_eigenvalues);
    }    

    return;
}

void WaveCurve::speed_map(std::vector<Curve> &arclength_speed, std::vector<Curve> &arclength_eigenvalues, std::vector<Curve> &arclength_reference_eigenvalues) const {
    arclength_speed.clear();
    arclength_eigenvalues.clear();
    arclength_reference_eigenvalues.clear();

    // First the second half-wavecurve will be processed. The arc length will be deemed negative, starting from the
    // initial point.
    //
    std::vector<Curve> temp_arclength_speed, temp_arclength_eigenvalues;

    half_speed_map(WAVECURVE_NEGATIVE_ARCLENGTH, beginnig_of_second_half, wavecurve.size(), temp_arclength_speed, temp_arclength_eigenvalues);

    // Revert and get the minimum arclength (it will be used later).
    //
    double min_arclength = temp_arclength_speed.back().curve.back()(0);

    for (int i = temp_arclength_speed.size() - 1; i >= 0; i--){
        arclength_speed.push_back(temp_arclength_speed[i]);
        arclength_eigenvalues.push_back(temp_arclength_eigenvalues[i]);
    }

    // Now compute the first half-wavecurve. Notice that half_speed_map will not clean the vectors
    // so it can be used like this. Also get the maximum arclength (it will be used later).
    //
    double max_arclength = arclength_speed.back().curve.back()(0);

    half_speed_map(WAVECURVE_POSITIVE_ARCLENGTH, 0, beginnig_of_second_half, arclength_speed, arclength_eigenvalues);

    // Add the reference eigenvalues (as curves).
    for (int i = 0; i < reference_point.e.size(); i++){
        Curve temp;
        RealVector p(2);

        p(1) = reference_point.e[i].r;
        p(0) = min_arclength;
        temp.curve.push_back(p);

        p(0) = arclength_speed.back().curve.back()(0);//max_arclength;
        temp.curve.push_back(p);

        arclength_reference_eigenvalues.push_back(temp);
    }

//    reference_eigenvalues.resize(reference_point.e.size());
//    for (int i = 0; i < reference_point.e.size(); i++) reference_eigenvalues(i) = reference_point.e[i].r;

    return;
}


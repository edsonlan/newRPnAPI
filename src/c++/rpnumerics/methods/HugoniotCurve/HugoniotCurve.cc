#include "HugoniotCurve.h"

HugoniotCurve::HugoniotCurve(const FluxFunction *ff, const AccumulationFunction *aa, const Boundary *bb) : f(ff), a(aa), boundary(bb), subphysics(0) {
    info_ = std::string("Abstract HugoniotCurve");

    method_ = UNDEFINED_HUGONIOT; // If this method shows up then the programmer forgot to set 
}

HugoniotCurve::HugoniotCurve(const FluxFunction *ff, const AccumulationFunction *aa, const Boundary *bb, SubPhysics *s) : f(ff), a(aa), boundary(bb), subphysics(s){
}

HugoniotCurve::~HugoniotCurve(){
}

void HugoniotCurve::curve(const ReferencePoint &ref, int type, std::vector<HugoniotPolyLine> &classified_curve){
    std::vector<Curve> cc;

    curve(ref, type, cc);

    // TEMPORARY
    for (int i = 0; i < cc.size(); i++){
        for (int j = 0; j < cc[i].curve.size(); j++){
            if (cc[i].curve[j].size() == 3) cc[i].curve[j](2) = 1.0;
        }
    }
    // TEMPORARY

    ColorCurve colorcurve(*f, *a);

    std::cout << "Will colorify." << std::endl;

    for (int i = 0; i < cc.size(); i++){

        std:cout << "Hugo., cc[i].curve is:" << std::endl;
        for (int j = 0; j < cc[i].curve.size(); j++) std::cout << cc[i].curve[j] << std::endl;

        if (cc[i].curve.size() > 0){

            std::vector<RealVector> segmented_curve;
            for (int j = 0; j < cc[i].curve.size() - 1; j++){
                segmented_curve.push_back(cc[i].curve[j]);
                segmented_curve.push_back(cc[i].curve[j + 1]);

                std::cout << "p will be " << cc[i].curve[j] << ", q will be " << cc[i].curve[j + 1] << std::endl;
            }

            std::vector<HugoniotPolyLine> temp_classified_curve;
            std::vector<RealVector> transition_list;

            std::cout << "    Classify begin, segmented_curve.size() = " << segmented_curve.size() << std::endl;
            colorcurve.classify_segmented_curve(segmented_curve, ref, temp_classified_curve, transition_list);
            std::cout << "    Classify end." << std::endl;

            for (int j = 0; j < temp_classified_curve.size(); j++) classified_curve.push_back(temp_classified_curve[j]);
        }
    }

    std::cout << "Colored ok." << std::endl;

    return;
}


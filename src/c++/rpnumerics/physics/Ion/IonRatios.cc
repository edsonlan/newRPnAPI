#include "IonRatios.h"

IonRatios::IonRatios(){
}

IonRatios::~IonRatios(){
}

// Only a function of v.
//
//int IonRatios::jet(const WaveState &w, JetMatrix &r, int degree) const {
//    if (degree >= 0){
//        double v = w(1);

//        r.set(0, 3.0);

//        if (degree >= 1){
//            r.set(0, 0, 0.0);

//            if (degree == 2){
//                r.set(0, 0, 0, 0.0);
//            }
//        }
//    }

//    return 2;
//}

int IonRatios::jet(const WaveState &w, JetMatrix &r, int degree) const {
    if (degree >= 0){
        double v = w(1);

        r.set(0, 3.0 + 0.5*v);

        if (degree >= 1){
            r.set(0, 0, 0.5);

            if (degree == 2){
                r.set(0, 0, 0, 0.0);
            }
        }
    }

    return 2;
}


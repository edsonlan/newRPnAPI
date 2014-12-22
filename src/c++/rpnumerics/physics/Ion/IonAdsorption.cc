#include "IonAdsorption.h"


IonAdsorption::IonAdsorption(){
}

IonAdsorption::~IonAdsorption(){
}

int IonAdsorption::jet(const WaveState &w, JetMatrix &a, int degree) const {
    if (degree >= 0){

        a.set(0, 0.0);

        if (degree >= 1){
            a.set(0, 0, 0.0);
            a.set(0, 1, 0.0);

            if (degree == 2){
                a.set(0, 0, 0, 0.0);
                a.set(0, 0, 1, 0.0);
                a.set(0, 1, 1, 0.0);
            }
        }
    }

    return 2;
}


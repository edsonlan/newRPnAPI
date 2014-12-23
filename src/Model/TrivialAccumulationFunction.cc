#include "TrivialAccumulationFunction.h"

TrivialAccumulationFunction::TrivialAccumulationFunction() : AccumulationFunction() {
}

TrivialAccumulationFunction::~TrivialAccumulationFunction(){
}

int TrivialAccumulationFunction::jet(const WaveState &w, JetMatrix &a, int degree) const {
    int n = w.stateSpaceDim();

    a.resize(n);

    if (degree >= 0){
        for (int i = 0; i < n; i++) a.set(i, w(i));

        if (degree >= 1){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < i; j++){
                    a.set(i, j, 0.0); 
                    a.set(j, i, 0.0);

                    a.set(i, i, 1.0);
                }
            }

            if (degree >= 2){
                for (int i = 0; i < n; i++){
                    for (int j = 0; j < n; j++){
                        for (int k = 0; k < n; k++) a.set(i, j, k, 0.0);
                    }
                }
            }
        }
    }


    return degree;
}

TrivialAccumulationFunction* TrivialAccumulationFunction::clone() const {
    return new TrivialAccumulationFunction;
}


#include "JDAccumulationFunction.h"

JDAccumulationFunction::JDAccumulationFunction() : AccumulationFunction() {
}

JDAccumulationFunction::~JDAccumulationFunction(){
}

int JDAccumulationFunction::jet(const WaveState &w, JetMatrix &f, int degree) const {
    f.resize(2);

    if (degree >= 0){
        double u = w(0);
        double v = w(1);

        f.set(0, u);
        f.set(1, v*u);

        if (degree >= 1){
            f.set(0, 0, 1.0);
            f.set(0, 1, 0.0);

            f.set(1, 0, v);
            f.set(1, 1, u);

            if (degree >= 2){
                // TODO: Check these ones.
                //
                f.set(0, 0, 0, 0.0);
                f.set(0, 0, 1, 0.0);
                f.set(0, 1, 0, 0.0);
                f.set(0, 1, 1, 0.0);

                f.set(1, 0, 0, 0.0);
                f.set(1, 0, 1, 1.0);
                f.set(1, 1, 0, 1.0);
                f.set(1, 1, 1, 0.0);
            }
        }
    }

    return degree;
}

JDAccumulationFunction* JDAccumulationFunction::clone() const {
    return new JDAccumulationFunction;
}


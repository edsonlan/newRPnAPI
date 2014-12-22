#include "JDFluxFunction.h"

JDFluxFunction::JDFluxFunction(Parameter *e) : FluxFunction() {
    epsilon_parameter = e;
}

JDFluxFunction::~JDFluxFunction(){
}

int JDFluxFunction::jet(const WaveState &w, JetMatrix &f, int degree) const {
    f.resize(2);

    double epsilon = epsilon_parameter->value();

    if (degree >= 0){
        double u = w(0);
        double v = w(1);

        f.set(0, (u - 1.0)*(u - 1.0) + (v - 1.0)*(v - 1.0));
        f.set(1, v*f.get(0) + epsilon*(v + 1.0)*(v + 1.0));

        if (degree >= 1){
            f.set(0, 0, 2.0*(u - 1.0));
            f.set(0, 1, 2.0*(v - 1.0));

            f.set(1, 0, 2.0*v*(u - 1.0));
            f.set(1, 1, (u - 1.0)*(u - 1.0) + (v - 1.0)*(v - 1.0) + 2.0*v*(v - 1.0) + epsilon*2.*(v + 1.0));

            if (degree >= 2){
                f.set(0, 0, 0, 2.0);
                f.set(0, 0, 1, 0.0);
                f.set(0, 1, 0, 0.0);
                f.set(0, 1, 1, 2.0);

                f.set(1, 0, 0, 2.0*v);
                f.set(1, 0, 1, 2.0*(u - 1.0));
                f.set(1, 1, 0, 2.0*(u - 1.0));
                f.set(1, 1, 1, 2.0*(v - 1.0) + 2.0*(v + (v - 1.0)) + 2.0*epsilon);
            }
        }
    }

    return degree;
}

double JDFluxFunction::alpha_dot(const RealVector &p) const {
    double v = p(1);
    double epsilon = epsilon_parameter->value();

    return 2.0*epsilon*(v + 1.0);
}

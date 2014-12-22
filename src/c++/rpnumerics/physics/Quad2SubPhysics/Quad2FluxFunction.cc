#include "Quad2FluxFunction.h"

Quad2FluxFunction::Quad2FluxFunction(Parameter *a1, Parameter *b1, Parameter *c1, Parameter *d1, Parameter *e1,
                                     Parameter *a2, Parameter *b2, Parameter *c2, Parameter *d2, Parameter *e2): FluxFunction(){
    a1_parameter = a1;
    b1_parameter = b1;
    c1_parameter = c1;
    d1_parameter = d1;
    e1_parameter = e1;

    a2_parameter = a2;
    b2_parameter = b2;
    c2_parameter = c2;
    d2_parameter = d2;
    e2_parameter = e2;
}


Quad2FluxFunction::~Quad2FluxFunction(){
}

int Quad2FluxFunction::jet(const WaveState & x, JetMatrix & y, int degree) const {
    y.resize(2);

    if (degree >= 0) {
        double u = x(0);
        double v = x(1);

        double a1 = a1_parameter->value();
        double b1 = b1_parameter->value();
        double c1 = c1_parameter->value();
        double d1 = d1_parameter->value();
        double e1 = e1_parameter->value();

        double a2 = a2_parameter->value();
        double b2 = b2_parameter->value();
        double c2 = c2_parameter->value();
        double d2 = d2_parameter->value();
        double e2 = e2_parameter->value();

        y.set(0, 0.5 * ((a1 * u * u) + (2.0 * b1 * u * v) + (c1 * v * v)) + d1 * u + e1*v);

        y.set(1, 0.5 * ((a2 * u * u) + (2.0 * b2 * u * v) + (c2 * v * v)) + d2 * u + e2*v);

        if (degree >= 1){
            y.set(0, 0, a1 * u + b1 * v + d1);
            y.set(0, 1, b1 * u + c1 * v + e1);

            y.set(1, 0, a2 * u + b2 * v + d2);
            y.set(1, 1, b2 * u + c2 * v + e2);

            if (degree >= 2) {
                y.set(0, 0, 0, a1);
                y.set(0, 1, 0, b1);
                y.set(0, 0, 1, b1);
                y.set(0, 1, 1, c1);

                y.set(1, 0, 0, a2);
                y.set(1, 1, 0, b2);
                y.set(1, 0, 1, b2);
                y.set(1, 1, 1, c2);
            }
        }
    }

    return 2;
}


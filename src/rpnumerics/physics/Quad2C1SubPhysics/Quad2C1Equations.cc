#include "Quad2C1Equations.h"

Quad2C1Equations::Quad2C1Equations() : Equations(){
    number_variables   = 3;
    number_equations   = 2;
    number_constraints = 1;

//    a0 = params.component(0);
//    b0 = params.component(1);
//    c0 = params.component(2);
//    d0 = params.component(3);
//    e0 = params.component(4);
//    f0 = params.component(5);
//    g0 = params.component(6);
//    h0 = params.component(7);
//    i0 = params.component(8);

//    a1 = params.component(9);
//    b1 = params.component(10);
//    c1 = params.component(11);
//    d1 = params.component(12);
//    e1 = params.component(13);
//    f1 = params.component(14);
//    g1 = params.component(15);
//    h1 = params.component(16);
//    i1 = params.component(17);

//    a2 = params.component(18);
//    b2 = params.component(19);
//    c2 = params.component(20);
//    d2 = params.component(21);
//    e2 = params.component(22);
//    f2 = params.component(23);
//    g2 = params.component(24);
//    h2 = params.component(25);
//    i2 = params.component(26);

    double a = -1.0;
    double b = 0.0;
    double c = 0.1;

    a0 = a;
    b0 = b;
    c0 = 1.0;
    d0 = 0.0;
    e0 = c;
        f0 = g0 = h0 = i0 = 0.0;

    a1 = b;
    b1 = 1.0;
    c1 = 0.0;
    d1 = -c;
    e1 = 0.0;
        f1 = g1 = h1 = i1 = 0.0;

    // f2 cannot be zero! Otherwise the Bhaskara won't work!
    //
    f2 = -1.0;
    a2 = b2 = c2 = g2 = h2 = 0.0;
    d2 = e2 = i2 = 1.0;
    j2 = 0.0;
}

Quad2C1Equations::~Quad2C1Equations(){
}

int Quad2C1Equations::obtain_W_from_U(const RealVector &p, double &V){
    double u = p(0);
    double v = p(1);

    double A = 0.5*f2;
    double B = g2*u + h2*v + i2;
    double C = 0.5*(a2*u*u +2*b2*u*v + c2*v*v) + d2*u + e2*v + j2;

    // Eq. A*w^2 + B*w + C = 0 is equivalent to out2 = 0 below.

    double w1, w2;
    int info = Utilities::Bhaskara(B/A, C/A, w1, w2);

    if (info != BHASKARA_COMPLEX_ROOTS){
        double w = p(2);

        if (std::abs(w - w1) < std::abs(w - w2)) V = w1;
        else                                     V = w2;

        return CONSTRAINED_VARIABLES_OK;
    }
    else return CONSTRAINED_VARIABLES_ERROR;
}

int Quad2C1Equations::compute(const RealVector &p, int degree, JetMatrix &Fjet, JetMatrix &Gjet, JetMatrix &Cjet){
    Fjet.resize(number_variables, number_equations);
    Gjet.resize(number_variables, number_equations);
    Cjet.resize(number_variables, number_constraints);

    double u = p(0);
    double v = p(1);

    double w;

    int info = obtain_W_from_U(p, w);
    if (info == CONSTRAINED_VARIABLES_ERROR) return COMPUTE_ERROR;

    if (degree >= 0){
        // TODO: What should be done if all the coefficients in out2 are zero??? This seems to imply that w is free.

        double out0 = 0.5*(a0*u*u +2.0*b0*u*v + c0*v*v) + d0*u + e0*v + 0.5*(f0*w*w + 2.0*g0*u*w +2.0*h0*v*w) + i0*w;
        double out1 = 0.5*(a1*u*u +2.0*b1*u*v + c1*v*v) + d1*u + e1*v + 0.5*(f1*w*w + 2.0*g1*u*w +2.0*h1*v*w) + i1*w;
        double out2 = 0.5*(a2*u*u +2.0*b2*u*v + c2*v*v) + d2*u + e2*v + 0.5*(f2*w*w + 2.0*g2*u*w +2.0*h2*v*w) + i2*w + j2;

        // Compute and fill F, G and C.
        //
        Fjet.set(0, out0);
        Fjet.set(1, out1);

        Gjet.set(0, u);
        Gjet.set(1, v);

        Cjet.set(0, out2); // By hypothesis.

        if (degree >= 1){
            double out00, out01, out02;
            double out10, out11, out12;
            double out20, out21, out22;

            out00 = a0*u + b0*v + d0 + g0*w;
            out01 = b0*u + c0*v + e0 + h0*w;
            out02 = f0*w + g0*u + h0*v + i0;

            out10 = a1*u + b1*v + d1 + g1*w;
            out11 = b1*u + c1*v + e1 + h1*w;
            out12 = f1*w + g1*u + h1*v + i1;

            out20 = a2*u + b2*v + d2 + g2*w;
            out21 = b2*u + c2*v + e2 + h2*w;
            out22 = f2*w + g2*u + h2*v + i2;

            // Flux.
            //
            Fjet.set(0, 0, out00);
            Fjet.set(0, 1, out01);
            Fjet.set(0, 2, out02);

            Fjet.set(1, 0, out10);
            Fjet.set(1, 1, out11);
            Fjet.set(1, 2, out12);

            // TODO: Trivial cases like this should be handled automagically.
            //
            Gjet.set(0, 0, 1.0);
            Gjet.set(0, 1, 0.0);
            Gjet.set(0, 2, 0.0);

            Gjet.set(1, 0, 0.0);
            Gjet.set(1, 1, 1.0);
            Gjet.set(1, 2, 0.0);

            // Constraints.
            //
            Cjet.set(0, 0, out20);
            Cjet.set(0, 1, out21);
            Cjet.set(0, 2, out22);

            if (degree >= 2){
                // Flux.
                //
                Fjet.set(0, 0, 0, a0);
                Fjet.set(0, 0, 1, b0);
                Fjet.set(0, 0, 2, g0);

                Fjet.set(0, 1, 0, b0);
                Fjet.set(0, 1, 1, c0);
                Fjet.set(0, 1, 2, h0);

                Fjet.set(0, 2, 0, g0);
                Fjet.set(0, 2, 1, h0);
                Fjet.set(0, 2, 2, f0);

                Fjet.set(1, 0, 0, a1);
                Fjet.set(1, 0, 1, b1);
                Fjet.set(1, 0, 2, g1);

                Fjet.set(1, 1, 0, b1);
                Fjet.set(1, 1, 1, c1);
                Fjet.set(1, 1, 2, h1);

                Fjet.set(1, 2, 0, g1);
                Fjet.set(1, 2, 1, h1);
                Fjet.set(1, 2, 2, f1);

                // Accumulation.
                //
                Gjet.set(0, 0, 0, 0.0);
                Gjet.set(0, 0, 1, 0.0);
                Gjet.set(0, 0, 2, 0.0);

                Gjet.set(0, 1, 0, 0.0);
                Gjet.set(0, 1, 1, 0.0);
                Gjet.set(0, 1, 2, 0.0);

                Gjet.set(0, 2, 0, 0.0);
                Gjet.set(0, 2, 1, 0.0);
                Gjet.set(0, 2, 2, 0.0);

                Gjet.set(1, 0, 0, 0.0);
                Gjet.set(1, 0, 1, 0.0);
                Gjet.set(1, 0, 2, 0.0);

                Gjet.set(1, 1, 0, 0.0);
                Gjet.set(1, 1, 1, 0.0);
                Gjet.set(1, 1, 2, 0.0);

                Gjet.set(1, 2, 0, 0.0);
                Gjet.set(1, 2, 1, 0.0);
                Gjet.set(1, 2, 2, 0.0);               

                // Constraints.
                //
                Cjet.set(0, 0, 0, a1);
                Cjet.set(0, 0, 1, b1);
                Cjet.set(0, 0, 2, g1);

                Cjet.set(0, 1, 0, b1);
                Cjet.set(0, 1, 1, c1);
                Cjet.set(0, 1, 2, h1);

                Cjet.set(0, 2, 0, g1);
                Cjet.set(0, 2, 1, h1);
                Cjet.set(0, 2, 2, f1);
            }
         }
    }

    return COMPUTE_OK;
}


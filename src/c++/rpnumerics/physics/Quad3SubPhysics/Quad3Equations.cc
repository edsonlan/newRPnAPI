#include "Quad3Equations.h"

Quad3Equations::Quad3Equations() : Equations(){
    number_variables   = 3;
    number_equations   = 3;
    number_constraints = 0;

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

    a2 = b2 = c2 = f2 = g2 = h2 = 0.0;
    d2 = e2 = i2 = 1.0;
    j2 = 0.0;
}

Quad3Equations::~Quad3Equations(){
}

int Quad3Equations::compute(const RealVector &p, int degree, JetMatrix &Fjet, JetMatrix &Gjet, JetMatrix &Cjet){
    Fjet.resize(number_variables, number_equations);
    Gjet.resize(number_variables, number_equations);
    Cjet.resize(0, 0);

    double u = p(0);
    double v = p(1);
    double w = p(2);

    if (degree >= 0){
        double out0 = 0.5*(a0*u*u +2*b0*u*v + c0*v*v) + d0*u + e0*v + 0.5*(f0*w*w + 2*g0*u*w +2*h0*v*w) + i0*w;
        double out1 = 0.5*(a1*u*u +2*b1*u*v + c1*v*v) + d1*u + e1*v + 0.5*(f1*w*w + 2*g1*u*w +2*h1*v*w) + i1*w;
        double out2 = 0.5*(a2*u*u +2*b2*u*v + c2*v*v) + d2*u + e2*v + 0.5*(f2*w*w + 2*g2*u*w +2*h2*v*w) + i2*w + j2;

        // Compute and fill F, G and C.
        //
        Fjet.set(0, out0);
        Fjet.set(1, out1);

        Gjet.set(0, u);
        Gjet.set(1, v);

        Fjet.set(2, out2);

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
            Fjet.set(2, 0, out20);
            Fjet.set(2, 1, out21);
            Fjet.set(2, 2, out22);

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
                Fjet.set(2, 0, 0, a1);
                Fjet.set(2, 0, 1, b1);
                Fjet.set(2, 0, 2, g1);

                Fjet.set(2, 1, 0, b1);
                Fjet.set(2, 1, 1, c1);
                Fjet.set(2, 1, 2, h1);

                Fjet.set(2, 2, 0, g1);
                Fjet.set(2, 2, 1, h1);
                Fjet.set(2, 2, 2, f1);
            }
         }
    }

    return COMPUTE_OK;
}



/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad3FluxFunction.cc
 **/

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "Quad3FluxFunction.h"
#include <math.h>

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

Quad3FluxFunction::Quad3FluxFunction(const Quad3FluxParams & params) : FluxFunction(params) {
}

Quad3FluxFunction::Quad3FluxFunction(const Quad3FluxFunction & copy) : FluxFunction(copy.fluxParams()) {
}

Quad3FluxFunction * Quad3FluxFunction::clone() const {
    return new Quad3FluxFunction(*this);
}

Quad3FluxFunction::~Quad3FluxFunction(void) {
}

// Quad3's jet function
//
// Given a vector x = (u, v, w) this function fills field, Jacobian and Hessian.
// For each dimension F is of the form:
//
//    0.5*(a*u^2 +2*b*u*v + c*v^2) + d*u + e*v + 0.5*(f*w^2 + 2*g*u*w +2*h*v*w) + i*w.
//
// Each row of the correspondent Jacobian is of the form:
//
//     [ dF/du    dF/dv    dF/dw ],
//
// where
//
//     dF/du = a*u + b*v + d + g*w,
//     dF/dv = b*u + c*v + e + h*w,
//     dF/dw = f*w + g*u + h*v + i.

int Quad3FluxFunction::jet(const WaveState & x, JetMatrix & y, int degree = 2) const {

    // Calculate F
    double a0, b0, c0, d0, e0, f0, g0, h0, i0;
    double a1, b1, c1, d1, e1, f1, g1, h1, i1;
    double a2, b2, c2, d2, e2, f2, g2, h2, i2;


    double out0, out1, out2;

    // Extract the parameters
    const FluxParams params = fluxParams();


    a0 = params.component(0);
    b0 = params.component(1);
    c0 = params.component(2);
    d0 = params.component(3);
    e0 = params.component(4);
    f0 = params.component(5);
    g0 = params.component(6);
    h0 = params.component(7);
    i0 = params.component(8);

    a1 = params.component(9);
    b1 = params.component(10);
    c1 = params.component(11);
    d1 = params.component(12);
    e1 = params.component(13);
    f1 = params.component(14);
    g1 = params.component(15);
    h1 = params.component(16);
    i1 = params.component(17);

    a2 = params.component(18);
    b2 = params.component(19);
    c2 = params.component(20);
    d2 = params.component(21);
    e2 = params.component(22);
    f2 = params.component(23);
    g2 = params.component(24);
    h2 = params.component(25);
    i2 = params.component(26);

    double u = x(0);
    double v = x(1);
    double w = x(2);

    // Compute and fill F
    out0 = 0.5*(a0*u*u +2*b0*u*v + c0*v*v) + d0*u + e0*v + 0.5*(f0*w*w + 2*g0*u*w +2*h0*v*w) + i0*w;
    out1 = 0.5*(a1*u*u +2*b1*u*v + c1*v*v) + d1*u + e1*v + 0.5*(f1*w*w + 2*g1*u*w +2*h1*v*w) + i1*w;
    out2 = 0.5*(a2*u*u +2*b2*u*v + c2*v*v) + d2*u + e2*v + 0.5*(f2*w*w + 2*g2*u*w +2*h2*v*w) + i2*w;

//    printf("Quad3 @ (%f, %f, %f) = (%f, %f, %f)\n", u, v, w, out0, out1, out2);

    y.set(0, out0);
    y.set(1, out1);
    y.set(2, out2);

    if (degree > 0) {
        // Calculate DF
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

        y.set(0, 0, out00);
        y.set(0, 1, out01);
        y.set(0, 2, out02);

        y.set(1, 0, out10);
        y.set(1, 1, out11);
        y.set(1, 2, out12);

        y.set(2, 0, out20);
        y.set(2, 1, out21);
        y.set(2, 2, out22);

        if (degree > 1) {
            // Calculate D2F
            y.set(0, 0, 0, a0);
            y.set(0, 0, 1, b0);
            y.set(0, 0, 2, g0);

            y.set(0, 1, 0, b0);
            y.set(0, 1, 1, c0);
            y.set(0, 1, 2, h0);

            y.set(0, 2, 0, g0);
            y.set(0, 2, 1, h0);
            y.set(0, 2, 2, f0);

            /**/

            y.set(1, 0, 0, a1);
            y.set(1, 0, 1, b1);
            y.set(1, 0, 2, g1);

            y.set(1, 1, 0, b1);
            y.set(1, 1, 1, c1);
            y.set(1, 1, 2, h1);

            y.set(1, 2, 0, g1);
            y.set(1, 2, 1, h1);
            y.set(1, 2, 2, f1);

            /**/

            y.set(2, 0, 0, a1);
            y.set(2, 0, 1, b1);
            y.set(2, 0, 2, g1);

            y.set(2, 1, 0, b1);
            y.set(2, 1, 1, c1);
            y.set(2, 1, 2, h1);

            y.set(2, 2, 0, g1);
            y.set(2, 2, 1, h1);
            y.set(2, 2, 2, f1);

            if (degree > 2) {
                return 0; //UNSUCCESSFUL_PROCEDURE;

            }
        }
    }

    return 2; //SUCCESSFUL_PROCEDURE;
}


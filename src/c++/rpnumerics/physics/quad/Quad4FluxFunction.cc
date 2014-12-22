
/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad2FluxFunction.cc
 **/

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "Quad4FluxFunction.h"
#include <math.h>

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

Quad4FluxFunction::Quad4FluxFunction(const Quad4FluxParams & params) : FluxFunction(params) {
}

Quad4FluxFunction::Quad4FluxFunction(const Quad4FluxFunction & copy) : FluxFunction(copy.fluxParams()) {
}

Quad4FluxFunction * Quad4FluxFunction::clone() const {
    return new Quad4FluxFunction(*this);
}

Quad4FluxFunction::~Quad4FluxFunction(void) {
}

//int Quad4FluxFunction::jet(const WaveState & x, JetMatrix & output, int degree = 2) const {
//
//    // Calculate F
//    // F = A + B x + (1/2) x^T C x
//
//    int spaceDim = 4;
//    double q [spaceDim];
//    double res [spaceDim];
//    for (int i = 0; i < spaceDim; i++)
//        q[i] = x(i);
//    for (int i = 0; i < spaceDim; i++) {
//        res[i] = Quad4FluxParams::DEFAULT_A[i];
//        for (int j = 0; j < spaceDim; j++) {
//            res[i] = res[i] + Quad4FluxParams::DEFAULT_B[i] [j] * q[j];
//            for (int k = 0; k < spaceDim; k++)
//                res[i] = res[i] + 0.5 * Quad4FluxParams::DEFAULT_C[i][j][k] * q[j] * q[k];
//        }
//    }
//    for (int i = 0; i < spaceDim; i++)
//        output(i, res[i]);
//
//    if (degree > 0) {
//        // Calculate DF
//        // DF = B + C x
//
//        double q [spaceDim];
//        double res [spaceDim][spaceDim];
//        for (int i = 0; i < spaceDim; i++)
//            q[i] = x(i);
//        for (int i = 0; i < spaceDim; i++)
//            for (int j = 0; j < spaceDim; j++) {
//                res[i] [j] = Quad4FluxParams::DEFAULT_B[i] [j];
//                for (int k = 0; k < spaceDim; k++)
//                    res[i] [j] = res[i] [j] + Quad4FluxParams::DEFAULT_C[i][j][k] * q[k];
//            }
//        for (int i = 0; i < spaceDim; i++)
//            for (int j = 0; j < spaceDim; j++)
//                output(i, j, res[i] [j]);
//
//
//        if (degree > 1) {
//
//            // Calculate D2F
//            // DFF = C
//            for (int i = 0; i < spaceDim; i++)
//                for (int j = 0; j < spaceDim; j++)
//                    for (int k = 0; k < spaceDim; k++)
//                        output(i, j, k, Quad4FluxParams::DEFAULT_C[i][j][k]);
//
//
//            if (degree > 2) {
//                return 0; //UNSUCCESSFUL_PROCEDURE;
//
//            }
//
//        }
//    }
//    return 2; //SUCCESSFUL_PROCEDURE;
//}

int Quad4FluxFunction::jet(const WaveState & x, JetMatrix & y, int degree = 2) const {


    // Calculate F
    double a0, b0, c0, d0, e0, f0, g0, h0, i0;
    double a1, b1, c1, d1, e1, f1, g1, h1, i1;
    double a2, b2, c2, d2, e2, f2, g2, h2, i2;
    double a3, b3, c3, d3, e3, f3, g3, h3, i3;

    double j0, k0, l0, m0, n0;
    double j1, k1, l1, m1, n1;
    double j2, k2, l2, m2, n2;
    double j3, k3, l3, m3, n3;

    double out0, out1, out2, out3;

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

    j0 = params.component(9);
    k0 = params.component(10);
    l0 = params.component(11);
    m0 = params.component(12);
    n0 = params.component(13);




    a1 = params.component(14);
    b1 = params.component(15);
    c1 = params.component(16);
    d1 = params.component(17);
    e1 = params.component(18);
    f1 = params.component(19);
    g1 = params.component(20);
    h1 = params.component(21);
    i1 = params.component(22);

    j1 = params.component(23);
    k1 = params.component(24);
    l1 = params.component(25);
    m1 = params.component(26);
    n1 = params.component(27);




    a2 = params.component(28);
    b2 = params.component(29);
    c2 = params.component(30);
    d2 = params.component(31);
    e2 = params.component(32);
    f2 = params.component(33);
    g2 = params.component(34);
    h2 = params.component(35);
    i2 = params.component(36);

    j2 = params.component(37);
    k2 = params.component(38);
    l2 = params.component(39);
    m2 = params.component(40);
    n2 = params.component(41);


    a3 = params.component(42);
    b3 = params.component(43);
    c3 = params.component(44);
    d3 = params.component(45);
    e3 = params.component(46);
    f3 = params.component(47);
    g3 = params.component(48);
    h3 = params.component(49);
    i3 = params.component(50);

    j3 = params.component(51);
    k3 = params.component(52);
    l3 = params.component(53);
    m3 = params.component(54);
    n3 = params.component(55);




    double u = x(0);
    double v = x(1);
    double w = x(2);
    double z = x(3);

    // Compute and fill F
    out0 = 0.5 * (a0 * u * u + 2 * b0 * u * v + c0 * v * v) + d0 * u + e0 * v + 0.5 * (f0 * w * w + 2 * g0 * u * w + 2 * h0 * v * w) + i0 * w + 0.5 * (j0 * z * z + 2 * k0 * u * z + 2 * l0 * v * z + 2 * m0 * w * z) + n0*z;

    out1 = 0.5 * (a1 * u * u + 2 * b1 * u * v + c1 * v * v) + d1 * u + e1 * v + 0.5 * (f1 * w * w + 2 * g1 * u * w + 2 * h1 * v * w) + i1 * w + 0.5 * (j1 * z * z + 2 * k1 * u * z + 2 * l1 * v * z + 2 * m1 * w * z) + n1*z;
    out2 = 0.5 * (a2 * u * u + 2 * b2 * u * v + c2 * v * v) + d2 * u + e2 * v + 0.5 * (f2 * w * w + 2 * g2 * u * w + 2 * h2 * v * w) + i2 * w + 0.5 * (j2 * z * z + 2 * k2 * u * z + 2 * l2 * v * z + 2 * m2 * w * z) + n2*z;
    out3 = 0.5 * (a3 * u * u + 2 * b3 * u * v + c3 * v * v) + d3 * u + e3 * v + 0.5 * (f3 * w * w + 2 * g3 * u * w + 2 * h3 * v * w) + i3 * w + 0.5 * (j3 * z * z + 2 * k3 * u * z + 2 * l3 * v * z + 2 * m3 * w * z) + n3*z;



    //    out1 = 0.5 * (a1 * u * u + 2 * b1 * u * v + c1 * v * v) + d1 * u + e1 * v + 0.5 * (f1 * w * w + 2 * g1 * u * w + 2 * h1 * v * w) + i1*w;
    //    out2 = 0.5*(a2*u*u +2*b2*u*v + c2*v*v) + d2*u + e2*v + 0.5*(f2*w*w + 2*g2*u*w +2*h2*v*w) + i2*w;
    //    out3 = 0.5*(a3*u*u +2*b3*u*v + c3*v*v) + d3*u + e3*v + 0.5*(f3*w*w + 2*g3*u*w +2*h3*v*w) + i3*w;


    //    printf("Quad3 @ (%f, %f, %f) = (%f, %f, %f)\n", u, v, w, out0, out1, out2);

    y.set(0, out0);
    y.set(1, out1);
    y.set(2, out2);
    y.set(3, out3);

    if (degree > 0) {
        // Calculate DF
        double out00, out01, out02, out03;
        double out10, out11, out12, out13;
        double out20, out21, out22, out23;
        double out30, out31, out32, out33;


        out00 = a0 * u + b0 * v + d0 + g0 * w + h0*z;
        out01 = b0 * u + c0 * v + e0 + h0 * w + l0*z;
        out02 = f0 * w + g0 * u + h0 * v + i0 + m0*z;
        out03 = j0 * z + k0 * u + l0 * v + m0 * w + n0;

        out10 = a1 * u + b1 * v + d1 + g1 * w + k1*z;
        out11 = b1 * u + c1 * v + e1 + h1 * w + l1*z;
        out12 = f1 * w + g1 * u + h1 * v + i1 + m1*z;
        out13 = j1 * z + k1 * u + l1 * v + m1 * w + n1;


        out20 = a2 * u + b2 * v + d2 + g2 * w + k2*z;
        out21 = b2 * u + c2 * v + e2 + h2 * w + l2*z;
        out22 = f2 * w + g2 * u + h2 * v + i2 + m2*z;
        out23 = j2 * z + k2 * u + l2 * v + m2 * w + n2;

        out30 = a3 * u + b3 * v + d3 + g3 * w + k3*z;
        out31 = b3 * u + c3 * v + e3 + h3 * w + l3*z;
        out32 = f3 * w + g3 * u + h3 * v + i3 + m3*z;
        out33 = j3 * z + k3 * u + l3 * v + m3 * w + n3;



        y.set(0, 0, out00);
        y.set(0, 1, out01);
        y.set(0, 2, out02);
        y.set(0, 3, out03);

        y.set(1, 0, out10);
        y.set(1, 1, out11);
        y.set(1, 2, out12);
        y.set(1, 3, out13);

        y.set(2, 0, out20);
        y.set(2, 1, out21);
        y.set(2, 2, out22);
        y.set(2, 3, out23);

        y.set(3, 0, out30);
        y.set(3, 1, out31);
        y.set(3, 2, out32);
        y.set(3, 3, out33);


        if (degree > 1) {
            // Calculate D2F
            y.set(0, 0, 0, a0);
            y.set(0, 0, 1, b0);
            y.set(0, 0, 2, g0);
            y.set(0, 0, 3, k0);

            y.set(0, 1, 0, b0);
            y.set(0, 1, 1, c0);
            y.set(0, 1, 2, h0);
            y.set(0, 1, 3, l0);



            y.set(0, 2, 0, g0);
            y.set(0, 2, 1, h0);
            y.set(0, 2, 2, f0);
            y.set(0, 2, 3, m0);
            
            
            
            y.set(0, 3, 0, k0);
            y.set(0, 3, 1, l0);
            y.set(0, 3, 2, m0);
            y.set(0, 3, 3, j0);


            /**/

            y.set(1, 0, 0, a1);
            y.set(1, 0, 1, b1);
            y.set(1, 0, 2, g1);
            y.set(1, 0, 3, k1);

            y.set(1, 1, 0, b1);
            y.set(1, 1, 1, c1);
            y.set(1, 1, 2, h1);
            y.set(1, 1, 3, l1);


            y.set(1, 2, 0, g1);
            y.set(1, 2, 1, h1);
            y.set(1, 2, 2, f1);
            y.set(1, 2, 3, m1);


            y.set(1, 3, 0, k1);
            y.set(1, 3, 1, l1);
            y.set(1, 3, 2, m1);
            y.set(1, 3, 3, j1);




            /**/

            y.set(2, 0, 0, a2);
            y.set(2, 0, 1, b2);
            y.set(2, 0, 2, g2);
            y.set(2, 0, 3, k2);

            y.set(2, 1, 0, b2);
            y.set(2, 1, 1, c2);
            y.set(2, 1, 2, h2);
            y.set(2, 1, 3, l2);


            y.set(2, 2, 0, g2);
            y.set(2, 2, 1, h2);
            y.set(2, 2, 2, f2);
            y.set(2, 2, 3, m2);


            y.set(2, 3, 0, k2);
            y.set(2, 3, 1, l2);
            y.set(2, 3, 2, m2);
            y.set(2, 3, 3, j2);


            /**/


            y.set(3, 0, 0, a3);
            y.set(3, 0, 1, b3);
            y.set(3, 0, 2, g3);
            y.set(3, 0, 3, k3);

            y.set(3, 1, 0, b3);
            y.set(3, 1, 1, c3);
            y.set(3, 1, 2, h3);
            y.set(3, 1, 3, l3);


            y.set(3, 2, 0, g3);
            y.set(3, 2, 1, h3);
            y.set(3, 2, 2, f3);
            y.set(3, 2, 3, m3);


            y.set(3, 3, 0, k3);
            y.set(3, 3, 1, l3);
            y.set(3, 3, 2, m3);
            y.set(3, 3, 3, j3);






            if (degree > 2) {
                return 0; //UNSUCCESSFUL_PROCEDURE;

            }
        }
    }

    return 2; //SUCCESSFUL_PROCEDURE;
}


























//
//    // Calculate F
//    double a1, b1, c1, d1, e1, a2, b2, c2, d2, e2, out0, out1;
//
//
//    const FluxParams params = fluxParams();
//
//    RealVector parVector = params.params();
//    
//
//    cout << parVector << "\n";
//
//
//    a1 = params.component(0);
//    b1 = params.component(1);
//    c1 = params.component(2);
//    d1 = params.component(3);
//    e1 = params.component(4);
//
//
//    cout << "a1 0: " << a1 << "\n";
//    cout << "b1 1: " << b1 << "\n";
//    cout << "c1 2: " << c1 << "\n";
//    cout << "d1 3: " << d1 << "\n";
//    cout << "e1 4: " << e1 << "\n";
//
//
//
//    a2 = params.component(5);
//    b2 = params.component(6);
//    c2 = params.component(7);
//    d2 = params.component(8);
//    e2 = params.component(9);
//    
//
//    cout << "a2 5: " << a2 << "\n";
//    cout << "b2 6: " << b2 << "\n";
//    cout << "c2 7: " << c2 << "\n";
//    cout << "d2 8: " << d2 << "\n";
//    cout << "e2 9: " << e2 << "\n";
//
//    
//    
//
////    double u = x(0);
////    double v = x(1);
//
//    double u = 0.14;
//    double v = 0.7;
//
//    
//    
//    out0 = 0.5 * (a1 * pow(u, (double) 2) + 2.0 * b1 * u * v + c1 * pow(v, (double) 2)) + d1 * u + e1*v;
//    
//    
//    
//    out1 = 0.5 * (a2 * pow(u, (double) 2) + 2.0 * b2 * u * v + c2 * pow(v, (double) 2)) + d2 * u + e2*v;


//out3 = 0.5 * (a3*u*u + 2 * b3*u*v + c2*v*v) + d2*u + e2*v + 0.5 * (f2*w*w + 2 * g2*u*w + 2 * h2*v*w) + i2*w;
//out3 = 0.5 * (a3*u*u + 2 * b3*u*v + c2*v*v) + d2*u + e2*v + 0.5 * (f2*w*w + 2 * g2*u*w + 2 * h2*v*w) + i2*w;
//
//    y.set(0, out0);
//    y.set(1, out1);
//
//
//    cout << "f (C++): " << y.set(0) << " "<< y.set(1)<<"\n";
//
//
//    
//    
//
//    if (degree > 0) {
//        // Calculate DF
//        double out00, out01, out10, out11;
//
//        out00 = a1 * u + b1 * v + d1;
//        out01 = b1 * u + c1 * v + e1;
//        out10 = a2 * u + b2 * v + d2;
//        out11 = b2 * u + c2 * v + e2;
//
//        y.set(0, 0, out00);
//        y.set(0, 1, out01);
//        y.set(1, 0, out10);
//        y.set(1, 1, out11);
//
//        if (degree > 1) {
//            // Calculate D2F
//            y.set(0, 0, 0, a1);
//            y.set(1, 0, 0, b1);
//            y.set(0, 1, 0, b1);
//            y.set(1, 1, 0, c1);
//            y.set(0, 0, 1, a2);
//            y.set(1, 0, 1, b2);
//            y.set(0, 1, 1, b2);
//            y.set(1, 1, 1, c2);
//
//            if (degree > 2) {
//                return 0; //UNSUCCESSFUL_PROCEDURE;
//
//            }
//        }
//    }
//
//    return 2; //SUCCESSFUL_PROCEDURE;
//}


























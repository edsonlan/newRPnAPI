#include "SimpleFlux.h"

void SimpleFlux::init(int n, const double *c){
    coef.resize(n);
    for (int i = 0; i < n; i++) coef(i) = c[i];

    coef_der1.resize(n - 1);
    for (int i = 0; i < n - 1; i++) coef_der1(i) = coef(i + 1)*(double)(i + 1);

    coef_der2.resize(n - 2);
    for (int i = 0; i < n - 2; i++) coef_der2(i) = coef(i + 2)*(double)(i + 2)*(double)(i + 1);

    std::cout << "coef      = " << coef << std::endl;
    std::cout << "coef_der1 = " << coef_der1 << std::endl;
    std::cout << "coef_der2 = " << coef_der2 << std::endl;

    return;
}

SimpleFlux::SimpleFlux(const RealVector &c) : FluxFunction() {
    init(c.size(), c.components());
}

SimpleFlux::SimpleFlux(int n, const double *c) : FluxFunction() {
    init(n, c);
}

SimpleFlux::~SimpleFlux(){
}

SimpleFlux* SimpleFlux::clone() const{
    return new SimpleFlux(coef);
}

int SimpleFlux::jet(const WaveState &w, JetMatrix &m, int degree) const {
    m.resize(2);

    double u = w(0);
    double v = 1000.0*w(1); //w(1);

//    std::cout << "u = " << u << ", v = " << v << std::endl;

    int n = coef.size();

    if (degree >= 0){
        RealVector x(n);
        x(0) = 1.0;

        // f = coef0 + u*coef1 + ... + u^(n - 1)*coef_{n - 1}
        for (int i = 1; i < n; i++) x(i) = x(i - 1)*u;

        m.set(0, x*coef);
        m.set(1, v);

        if (degree >= 1){
//            RealVector coef_der1(n - 1);
//            for (int i = 0; i < n - 1; i++) coef_der1(i) = coef(i + 1)*(double)(i + 1);
//            std::cout << "For Jacobian: coef_der1 = " << coef_der1 << std::endl;
            
            // df_du = coef1 + ... + (n - 1)*u^(n - 2)*coef_{n - 1}
            double df_du = 0.0;
            for (int i = 0; i < n - 1; i++) df_du += coef_der1(i)*x(i);
//            df_du = coef(1) + 2.0*u*coef(2);

            m.set(0, 0, df_du); // dF0/du
            m.set(0, 1, 0.0);   // dF0/dv = 0.0
            m.set(1, 0, 0.0);   // dF1/du = 0.0
            m.set(1, 1, 1000.0);   // dF1/dv = 1.0

            if (degree >= 2){
//                RealVector coef_der2(n - 2);                
//                for (int i = 0; i < n - 2; i++) coef_der2(i) = coef_der1(i + 1)*(double)(i + 2);

                double d2f_du2 = 0.0;
                for (int i = 0; i < n - 2; i++) d2f_du2 += coef_der2(i)*x(i);

                m.set(0, 0, 0, d2f_du2);
                m.set(0, 0, 1, 0.0);
                m.set(0, 1, 0, 0.0);
                m.set(0, 1, 1, 0.0);

                m.set(1, 0, 0, 0.0);
                m.set(1, 0, 1, 0.0);
                m.set(1, 1, 0, 0.0);
                m.set(1, 1, 1, 0.0);
            }
        }
    }

    return 2;
}


#ifndef _QUAD3EQUATIONS_
#define _QUAD3EQUATIONS_

#include "Equations.h"

class Quad3Equations : public Equations {
    private:
    protected:
        double a0, b0, c0, d0, e0, f0, g0, h0, i0;
        double a1, b1, c1, d1, e1, f1, g1, h1, i1;
        double a2, b2, c2, d2, e2, f2, g2, h2, i2, j2;

        int compute(const RealVector &p, int degree, JetMatrix &Fjet, JetMatrix &Gjet, JetMatrix &Cjet);
    public:
        Quad3Equations();
        virtual ~Quad3Equations();
};

#endif // _QUAD2C1EQUATIONS_


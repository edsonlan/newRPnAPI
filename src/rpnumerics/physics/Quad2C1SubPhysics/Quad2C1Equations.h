#ifndef _QUAD2C1EQUATIONS_
#define _QUAD2C1EQUATIONS_

#include "Equations.h"
#include "Utilities.h"

class Quad2C1Equations : public Equations {
    private:
    protected:
        double a0, b0, c0, d0, e0, f0, g0, h0, i0;
        double a1, b1, c1, d1, e1, f1, g1, h1, i1;
        double a2, b2, c2, d2, e2, f2, g2, h2, i2, j2;

    public:
        Quad2C1Equations();
        virtual ~Quad2C1Equations();

        int obtain_W_from_U(const RealVector &p, double &W);
        int compute(const RealVector &p, int degree, JetMatrix &Fjet, JetMatrix &Gjet, JetMatrix &Cjet);
};

#endif // _QUAD2C1EQUATIONS_


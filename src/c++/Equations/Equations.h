#ifndef _EQUATIONS_
#define _EQUATIONS_

#define CONSTRAINEDCONVECTION_OK    0
#define CONSTRAINEDCONVECTION_ERROR 1

#define CONSTRAINED_VARIABLES_OK      0
#define CONSTRAINED_VARIABLES_ERROR (-1)

#define COMPUTE_OK      0
#define COMPUTE_ERROR (-1)

#include "RealVector.h"
#include "DoubleMatrix.h"
#include "JetMatrix.h"

class Equations {
    private:
    protected:
        virtual int compute(const RealVector &p, int degree, JetMatrix &Fjet, JetMatrix &Gjet, JetMatrix &Cjet) = 0;

        // Derived classes must set these.
        //
        int number_variables; // number_variables = number_equations + number_constraints.
        int number_equations;
        int number_constraints;

        virtual int obtain_W_from_U(const RealVector &p, RealVector &W){return CONSTRAINED_VARIABLES_OK;}
    public:
        Equations();
        virtual ~Equations();

        virtual int constrained_conservation_laws_jet(const RealVector &p, DoubleMatrix &DFcompleted, DoubleMatrix &DG);
};

#endif // _EQUATIONS_


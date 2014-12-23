#ifndef _ICDOWCHEMISTRY_
#define _ICDOWCHEMISTRY_

#include "JetMatrix.h"

class ICDOWChemistry {
    private:
    protected:
    public:
        ICDOWChemistry();
        virtual ~ICDOWChemistry();

        void aqueous_carbon_excess_oxygen(double hyd, int degree, JetMatrix &rho_CO_jet);

        void aqueous_hydrogen(double hyd, int degree, JetMatrix &rho_H_jet);

        void carbon_in_oil(double hyd, int degree, JetMatrix &rho_C_oil_jet);
};

#endif // _ICDOWCHEMISTRY_


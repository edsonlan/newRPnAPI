#ifndef _POLYDISPERSE_PARAMS_
#define _POLYDISPERSE_PARAMS_

#include "FluxParams.h"

class Polydisperse_Params : public FluxParams {
    private:
    protected:
    public:
        Polydisperse_Params(const double phimax,
                              const double rho1, const double rho2,
                              const double dVinf1, const double dVinf2,
                              const double n1, const double n2);
        Polydisperse_Params();

        Polydisperse_Params(const RealVector &);

        ~Polydisperse_Params();

};

#endif // _POLYDISPERSE_PARAMS_

#ifndef _COINCIDENCE_
#define _COINCIDENCE_

#include "RealVector.h"
#include <limits>

// Base class for the coincidence.
//
class Coincidence {
    private:
    protected:
    public:
        virtual ~Coincidence(){}

        virtual void lambdas(const RealVector &u, double &lambda_s, double &lambda_e, double &lambda_diff) const = 0;

        virtual double lambda_s(const RealVector &p) const = 0;
        virtual double lambda_e(const RealVector &p) const = 0;
        virtual double lambda_diff(const RealVector &p) const = 0;

        // For the Extension Map:
        //
        virtual bool extension_basis(const RealVector &u, double &fe, double &se) const = 0;
};

#endif // _COINCIDENCE_


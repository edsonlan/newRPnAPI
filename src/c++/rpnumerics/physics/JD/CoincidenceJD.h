#ifndef _COINCIDENCEJD_
#define _COINCIDENCEJD_

#include "Coincidence.h"
#include "JDFluxFunction.h"
#include "JDAccumulationFunction.h"

class CoincidenceJD : public Coincidence {
    private:
    protected:
        const JDFluxFunction         *flux;
        const JDAccumulationFunction *accum;
    public:
        CoincidenceJD(const JDFluxFunction *f, const JDAccumulationFunction *a);
        virtual ~CoincidenceJD();

        void lambdas(const RealVector &u, double &lambda_s, double &lambda_e, double &lambda_diff) const;

        double lambda_s(const RealVector &p) const;
        double lambda_e(const RealVector &p) const;
        double lambda_diff(const RealVector &p) const;

        // For the Extension Map:
        //
        bool extension_basis(const RealVector &u, double &fe, double &se) const;
};

#endif // _COINCIDENCEJD_


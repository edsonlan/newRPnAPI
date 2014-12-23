#ifndef _ICDOWCOINCIDENCE_
#define _ICDOWCOINCIDENCE_

#include "Coincidence.h"
#include "ICDOWHydrodynamics.h"

class ICDOWCoincidence : public Coincidence {
    private:
    protected:
        ICDOWHydrodynamics  *hydro;

        Parameter *phi_parameter;
    public:
        ICDOWCoincidence(ICDOWHydrodynamics *hy, Parameter *phi);
        virtual ~ICDOWCoincidence();

        void lambdas(const RealVector &p, double &lambda_s, double &lambda_e, double &lambda_diff) const;

        double lambda_s(const RealVector &p) const;
        double lambda_e(const RealVector &p) const;
        double lambda_diff(const RealVector &p) const;

        // For the Extension Map:
        //
        bool extension_basis(const RealVector &u, double &fe, double &se) const;
};

#endif // _ICDOW_


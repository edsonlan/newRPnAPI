#ifndef _ICDOWFLUXFUNCTION_
#define _ICDOWFLUXFUNCTION_

#include "FluxFunction.h"
#include "ICDOWChemistry.h"
#include "ICDOWHydrodynamics.h"

#define JETTESTER_ENABLED_ICDOWFLUX



class ICDOWFluxFunction: public FluxFunction
                         
{
    private:
    protected:
        ICDOWChemistry      *chemistry;
        ICDOWHydrodynamics  *hydro;
    public:
        ICDOWFluxFunction(ICDOWChemistry *ch, ICDOWHydrodynamics *hy);
        virtual ~ICDOWFluxFunction();

        int reduced_jet(const WaveState &u, JetMatrix &m, int degree) const;
        int reduced_jet(const RealVector &u, JetMatrix &m, int degree) const {
            WaveState w(u);

            int info = reduced_jet(w, m, degree);

            return info;
        }

        int jet(const WaveState &u, JetMatrix &m, int degree) const;

       
};

#endif // _ICDOWFLUXFUNCTION_


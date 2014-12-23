#ifndef _DEADVOLATILEVOLATILEGASCOINCIDENCE_
#define _DEADVOLATILEVOLATILEGASCOINCIDENCE_

#include "Coincidence.h"
#include "DeadVolatileVolatileGasThermodynamics.h"
#include "DeadVolatileVolatileGasHydrodynamics.h"

class DeadVolatileVolatileGasCoincidence : public Coincidence {
    private:
    protected:
        DeadVolatileVolatileGasThermodynamics *thermo;
        DeadVolatileVolatileGasHydrodynamics  *hydro;

        Parameter *re_parameter, *rg_parameter, *phi_parameter;
    public:
        DeadVolatileVolatileGasCoincidence(DeadVolatileVolatileGasThermodynamics *th, 
                                   DeadVolatileVolatileGasHydrodynamics *hy, 
                                   Parameter *re, Parameter *rg, 
                                   Parameter *phi);

        virtual ~DeadVolatileVolatileGasCoincidence();

        void lambdas(const RealVector &p, double &lambda_s, double &lambda_e, double &lambda_diff) const;

        double lambda_s(const RealVector &p) const;
        double lambda_e(const RealVector &p) const;
        double lambda_diff(const RealVector &p) const;

        // For the Extension Map:
        //
        bool extension_basis(const RealVector &u, double &fe, double &se) const;

        // To be used by DeadVolatileVolatileGasEvaporationExtension.
        //
        DeadVolatileVolatileGasThermodynamics* thermodynamics(){return thermo;}
};

#endif // _DEADVOLATILEVOLATILEGASCOINCIDENCE_


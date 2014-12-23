#ifndef _DEADVOLATILEVOLATILEGASEVAPORATIONEXTENSION_
#define _DEADVOLATILEVOLATILEGASEVAPORATIONEXTENSION_

#include "Extension.h"
#include "DeadVolatileVolatileGasCoincidence.h"
#include "DeadVolatileVolatileGasFluxFunction.h"
#include "DeadVolatileVolatileGasAccumulationFunction.h"
#include "Utilities.h"

class DeadVolatileVolatileGasEvaporationExtension : public Extension {
    private:
    protected:
        const DeadVolatileVolatileGasFluxFunction *flux;
        const DeadVolatileVolatileGasAccumulationFunction *accumulation;
        DeadVolatileVolatileGasCoincidence *god;
        Parameter *phi_parameter;
    public:
        DeadVolatileVolatileGasEvaporationExtension(const DeadVolatileVolatileGasFluxFunction *f, const DeadVolatileVolatileGasAccumulationFunction *a, DeadVolatileVolatileGasCoincidence *g, Parameter *phi);
        ~DeadVolatileVolatileGasEvaporationExtension();

        int extension(const RealVector &p, RealVector &ext_p);
        std::string name() const {return std::string("DeadVolatileVolatileGasEvaporationExtension");}
        int extension_type(){return EXPLICIT_EXTENSION;}
};

#endif // _DEADVOLATILEVOLATILEGASEVAPORATIONEXTENSION_


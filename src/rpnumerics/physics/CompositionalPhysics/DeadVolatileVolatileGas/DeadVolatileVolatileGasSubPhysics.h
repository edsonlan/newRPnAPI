#ifndef _DEADVOLATILEVOLATILEGASSUBPHYSICS_
#define _DEADVOLATILEVOLATILEGASSUBPHYSICS_

#include "SubPhysics.h"
#include "DeadVolatileVolatileGasThermodynamics.h"
#include "DeadVolatileVolatileGasHydrodynamics.h"
#include "DeadVolatileVolatileGasCoincidence.h"
#include "DeadVolatileVolatileGasEvaporationExtension.h"

#include "DeadVolatileVolatileGasFluxFunction.h"
#include "DeadVolatileVolatileGasAccumulationFunction.h"
#include "DeadVolatileVolatileGasCompositeCurve.h"
#include "LSODE.h"

#include "RectBoundary.h"

#include "Hugoniot_TP.h"
//#include "RarefactionCurve.h"

#include "HugoniotContinuation3D2D.h"
#include "HugoniotContinuation_nDnD.h"

#define DEADVOLATILEVOLATILEGAS_GENERIC_POINT 0

class DeadVolatileVolatileGasSubPhysics : public SubPhysics {
    private:
    protected:
        Parameter *B_parameter, *D_parameter, *mu_oB_parameter, *mu_oD_parameter, *mu_G_parameter;
        DeadVolatileVolatileGasThermodynamics *thermo;

        DeadVolatileVolatileGasHydrodynamics  *hydro;

        Parameter *re_parameter, *rg_parameter, *phi_parameter;

        DeadVolatileVolatileGasEvaporationExtension *evap_;
    public:
        DeadVolatileVolatileGasSubPhysics();
        virtual ~DeadVolatileVolatileGasSubPhysics();

        void shock_cases(std::vector<int> &type, std::vector<std::string> &name) const;
};

#endif // _DEADVOLATILEVOLATILEGASSUBPHYSICS_

#ifndef _TPCWSUBPHYSICS_
#define _TPCWSUBPHYSICS_

#include "SubPhysics.h"
#include "MolarDensity.h"
#include "VLE_Flash_TPCW.h"
#include "Thermodynamics.h"

#include "Accum2Comp2PhasesAdimensionalized.h"
#include "Flux2Comp2PhasesAdimensionalized.h"
#include "RectBoundary.h"
#include "Hugoniot_TP.h"

class TPCWSubPhysics: public SubPhysics {
    private:
    protected:
        // Flux parameters.
        //
        Parameter *abs_perm_parameter, *sin_beta_parameter;
        Parameter *cnw_parameter,  *cng_parameter;
        Parameter *expw_parameter, *expg_parameter;

        // Accumulation parameter.
        //
        Parameter *phi_parameter;

        MolarDensity *mdv, *mdl;
        VLE_Flash_TPCW *flash;
        Thermodynamics *tc;
    public:
        TPCWSubPhysics();
        virtual ~TPCWSubPhysics();
};

#endif // _TPCWSUBPHYSICS_


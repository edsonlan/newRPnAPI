#ifndef _HYSTERESIS_
#define _HYSTERESIS_

#include "Extension_Curve.h"
#include "Inflection_Curve.h"
#include "Boundary.h"

// This class implements the Hysteresis, which is the extension 
// of the inflection curve.
//
class Hysteresis {
    private:
    protected:
    public:
    static void curve(const FluxFunction *curve_flux,
            const AccumulationFunction *curve_accum,GridValues &,
            int characteristic_where, int curve_family,
            int domain_family,
            const FluxFunction *domain_ff,
            const AccumulationFunction *domain_aa,
            int singular,
            std::vector<RealVector> &curve_segments,
            std::vector<RealVector> &domain_segments);
};

#endif // _HYSTERESIS_


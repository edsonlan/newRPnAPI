#include "Hysteresis.h"

void Hysteresis::curve(
        const FluxFunction *curve_flux, const AccumulationFunction *curve_accum, GridValues &gv,
        int characteristic_where, int curve_family,
        int domain_family,
        const FluxFunction *domain_ff, const AccumulationFunction *domain_aa,
        int singular,
        std::vector<RealVector> &curve_segments,
        std::vector<RealVector> &domain_segments) {

    // Inflection curve

    std::vector<RealVector> inflectionSegments;

    Inflection_Curve inflectionCurve;
    inflectionCurve.curve(curve_flux, curve_accum, gv, curve_family, inflectionSegments);

    if (inflectionSegments.size() < 2) return;

    // Verify that the points in the inflection curve have the same dimension as the grid.
    if (inflectionSegments[0].size() < gv.grid(0).size()){
        int n = gv.grid(0).size();
        int m = inflectionSegments[0].size();

        for (int i = 0; i < inflectionSegments.size(); i++){
            inflectionSegments[i].resize(n);
            for (int j = m; j < n; j++) inflectionSegments[i].component(j) = 1.0;
        }
    }

    // Compute the extension curve for the inflection curve.
    Extension_Curve extension_curve;

    extension_curve.curve(domain_ff, domain_aa, curve_flux, curve_accum, gv, characteristic_where, singular, curve_family,
                          inflectionSegments, curve_segments, domain_segments);


    return;
}


/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) EllipticExtension.h
 */

#ifndef _EllipticExtension_H
#define _EllipticExtension_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "Elliptic_Boundary.h"
#include "Extension_Curve.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class EllipticExtension {
private:

public:

    int curve(const FluxFunction *f, const AccumulationFunction *a, int where_is_characteristic, int family,
            GridValues &g, std::vector<RealVector> &elliptic_extension_on_curve, std::vector<RealVector> &elliptic_extension_on_domain);

};

#endif //! _EllipticExtension_H

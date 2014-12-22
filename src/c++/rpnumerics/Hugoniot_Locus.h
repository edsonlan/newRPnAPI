/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Hugoniot_Locus.h
 */

#ifndef _Hugoniot_Locus_H
#define _Hugoniot_Locus_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */

#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Matrix.h"
#include "RealVector.h"

#include "ImplicitFunction.h"
#include "ContourMethod.h"
#include "ColorCurve.h"

#include "ReferencePoint.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class Hugoniot_Locus : public ImplicitFunction {
private:

 

public:


//    virtual int curve(const FluxFunction *f, const AccumulationFunction *a,
//            GridValues &g, const RealVector &r,
//            std::vector<RealVector> &hugoniot_curve)=0;


    virtual int classified_curve(const FluxFunction *f, const AccumulationFunction *a,
            GridValues &g, const ReferencePoint &r,
            std::vector<HugoniotPolyLine> &hugoniot_curve) = 0;

    virtual int classified_curve(const FluxFunction *f, const AccumulationFunction *a, 
                             GridValues &g, const ReferencePoint &r, 
                             std::vector<HugoniotPolyLine> &hugoniot_curve,
                             std::vector<bool> &circular)=0;

};

#endif //! _Hugoniot_Locus_H

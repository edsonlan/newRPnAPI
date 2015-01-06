#ifndef _REFERENCEPOINT_
#define _REFERENCEPOINT_

#include <stdio.h>
#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Viscosity_Matrix.h"
#include "RealVector.h"
#include "eigen.h"

#include "Matrix.h"
#include "DoubleMatrix.h"

class ReferencePoint {
    private:
        const FluxFunction *f;
        const AccumulationFunction *g;
        const Viscosity_Matrix *v;

    protected:
    public:
        //
        RealVector point;

        RealVector F, G;
        DoubleMatrix JF, JG;
        ViscosityJetMatrix B;

        std::vector< eigenpair > e;

        // Default constructor.
        ReferencePoint();
        
        ReferencePoint (const ReferencePoint &);

        // Constructor. Fills an inital point.
        //
        ReferencePoint(const RealVector &p,
                       const FluxFunction *ff, const AccumulationFunction *gg,
                       const Viscosity_Matrix *vv);

        // Destructor. Does nothing so far.
        //
        ~ReferencePoint(){}

        void fill_point(const RealVector &p,
                        const FluxFunction *ff, const AccumulationFunction *gg,
                        const Viscosity_Matrix *vv);

        void fill_point(const RealVector &p);

        ReferencePoint operator=(const ReferencePoint &orig);
};

#endif // _REFERENCEPOINT_


#ifndef _QUAD2SUBPHYSICS_
#define _QUAD2SUBPHYSICS_

#include "SubPhysics.h"
#include "Quad2FluxFunction.h"
#include "Quad2AccumulationFunction.h"
#include "LSODE.h"
#include "ImplicitHugoniotCurve.h"
#include "HugoniotContinuation2D2D.h"
#include "WaveCurveFactory.h"
#include "RectBoundary.h"
#include "Implicit_Extension_Curve.h"
#include "Quad2ExplicitHugoniotCurve.h"

class Quad2SubPhysics: public SubPhysics {
    private:
    protected:
        Parameter *a1_parameter, *b1_parameter, *c1_parameter, *d1_parameter, *e1_parameter;
        Parameter *a2_parameter, *b2_parameter, *c2_parameter, *d2_parameter, *e2_parameter;

    public:
        Quad2SubPhysics();
        virtual ~Quad2SubPhysics();

        double a1(){return a1_parameter->value();}
        double b1(){return b1_parameter->value();}
        double c1(){return c1_parameter->value();}
        double d1(){return d1_parameter->value();}
        double e1(){return e1_parameter->value();}

        double a2(){return a2_parameter->value();}
        double b2(){return b2_parameter->value();}
        double c2(){return c2_parameter->value();}
        double d2(){return d2_parameter->value();}
        double e2(){return e2_parameter->value();}
};

#endif // _QUAD2SUBPHYSICS_


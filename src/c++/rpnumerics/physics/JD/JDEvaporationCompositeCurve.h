#ifndef _JDEVAPORATIONCOMPOSITECURVE_
#define _JDEVAPORATIONCOMPOSITECURVE_

#include "CompositeCurve.h"
#include "JDEvap_Extension.h"

class JDEvaporationCompositeCurve : public CompositeCurve {
    private:
    protected:
        JDEvap_Extension *evap;
    public:
        JDEvaporationCompositeCurve(const AccumulationFunction *a, const FluxFunction *f, const Boundary *b, JDEvap_Extension *e);
        virtual ~JDEvaporationCompositeCurve();

        virtual int curve(const AccumulationFunction *RarAccum, const FluxFunction *RarFlux,
                          const Boundary *RarBoundary, 
                          const Curve &rarcurve,
                          const RealVector &composite_initial_point,
                          int last_point_in_rarefaction,
                          const ODE_Solver *odesolver,
                          double deltaxi,
                          int where_composite_begins, int fam, 
                          Curve &new_rarcurve,
                          Curve &compositecurve,
                          RealVector &final_direction,
                          int &reason_why,
                          int &edge);

};

#endif // _JDEVAPORATIONCOMPOSITECURVE_


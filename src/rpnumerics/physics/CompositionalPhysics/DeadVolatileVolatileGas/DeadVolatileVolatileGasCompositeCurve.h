#ifndef _DEADVOLATILEVOLATILEGASCOMPOSITECURVE_
#define _DEADVOLATILEVOLATILEGASCOMPOSITECURVE_

#include "CompositeCurve.h"
#include "DeadVolatileVolatileGasEvaporationExtension.h"

class DeadVolatileVolatileGasCompositeCurve: public CompositeCurve {
    private:
    protected:
        DeadVolatileVolatileGasEvaporationExtension *evapext;
    public:
        DeadVolatileVolatileGasCompositeCurve(DeadVolatileVolatileGasEvaporationExtension *e, 
                                      const AccumulationFunction *a, 
                                      const FluxFunction *f, 
                                      const Boundary *b, 
                                      ShockCurve *s, 
                                      Explicit_Bifurcation_Curves *ebc);

        virtual ~DeadVolatileVolatileGasCompositeCurve();

        int curve(const AccumulationFunction *RarAccum, const FluxFunction *RarFlux,
                  const Boundary *RarBoundary, 
                  const Curve &rarcurve,
                  const RealVector &composite_initial_point,
                  int last_point_in_rarefaction,
                  const ODE_Solver *odesolver,
                  double deltaxi,
                  void *linobj, double (*linear_function)(void *o, const RealVector &p),
                  int where_composite_begins, int fam, 
                  Curve &new_rarcurve,
                  Curve &compositecurve,
                  RealVector &final_direction,
                  int &reason_why,
                  int &edge);
};

#endif // _DEADVOLATILEVOLATILEGASCOMPOSITECURVE_


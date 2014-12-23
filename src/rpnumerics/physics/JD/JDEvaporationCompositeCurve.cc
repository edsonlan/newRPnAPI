#include "JDEvaporationCompositeCurve.h"

JDEvaporationCompositeCurve::JDEvaporationCompositeCurve(const AccumulationFunction *a, 
                                                         const FluxFunction *f, 
                                                         const Boundary *b, 
                                                         JDEvap_Extension *e) : CompositeCurve(a, f, b, 0, 0), evap(e){
}

JDEvaporationCompositeCurve::~JDEvaporationCompositeCurve(){
}

int JDEvaporationCompositeCurve::curve(const AccumulationFunction *RarAccum, const FluxFunction *RarFlux,
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
                                       int &edge){

    compositecurve.type = COMPOSITE_CURVE;
    compositecurve.family = fam;

    const Coincidence *coincidence = evap->coincidence();

    // TODO: Only works after the evaporation rarefaction reached the inflection and the coincidence.
    // 
    if (where_composite_begins == COMPOSITE_BEGINS_AT_INFLECTION){
        add_point_to_curve(rarcurve.curve.back(), rarcurve.curve.size() - 1, rarcurve, compositecurve);

        for (int i = rarcurve.curve.size() - 2; i >= 0; i--){
             // Only extensions from a evaporation curve are acceptable.
             //
             double lambda_s = coincidence->lambda_s(rarcurve.curve[i]);
             double lambda_e = coincidence->lambda_e(rarcurve.curve[i]);

             if (fam == 0 && (lambda_s < lambda_e)) return COMPOSITE_ERROR;
             if (fam == 1 && (lambda_s > lambda_e)) return COMPOSITE_ERROR;

             RealVector ext_p;
             int info_ext = evap->extension(rarcurve.curve[i], ext_p);

             // TODO: Check this.
             //
             if (info_ext == EXTENSION_ERROR) {
                 reason_why = EXTENSION_ERROR;
                 return COMPOSITE_ERROR;
             }

             if (!boundary->inside(ext_p)){
                 RealVector r;
                 boundary->intersection(compositecurve.curve.back(), ext_p, r, edge);

                 compositecurve.curve.push_back(r);

                 reason_why = COMPOSITE_REACHED_BOUNDARY;
                 return COMPOSITE_OK;
             }

             add_point_to_curve(ext_p, i, rarcurve, compositecurve);
        }

        reason_why = COMPOSITE_COMPLETED;
    }
    else {
        // TODO: COMPOSITE_AFTER_COMPOSITE, nothing to do so far.
    }
    
    if (compositecurve.curve.size() > 1){
        int n = compositecurve.curve.size();

        final_direction = compositecurve.curve[n - 1] - compositecurve.curve[n - 2];
        normalize(final_direction);

        compositecurve.final_direction = final_direction;
        compositecurve.last_point = compositecurve.curve.back();
        compositecurve.reason_to_stop = reason_why;
    }

    return COMPOSITE_OK;
}


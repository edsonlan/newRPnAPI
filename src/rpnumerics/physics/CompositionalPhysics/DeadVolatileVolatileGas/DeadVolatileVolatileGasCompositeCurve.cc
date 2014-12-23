#include "DeadVolatileVolatileGasCompositeCurve.h"

DeadVolatileVolatileGasCompositeCurve::DeadVolatileVolatileGasCompositeCurve(DeadVolatileVolatileGasEvaporationExtension *e,
                                                             const AccumulationFunction *a, 
                                                             const FluxFunction *f, 
                                                             const Boundary *b, 
                                                             ShockCurve *s, 
                                                             Explicit_Bifurcation_Curves *ebc): CompositeCurve(a, f, b, s, ebc), evapext(e){

}

DeadVolatileVolatileGasCompositeCurve::~DeadVolatileVolatileGasCompositeCurve(){
}

// TODO: Only deals with COMPOSITE_BEGINS_AT_INFLECTION.
//
int DeadVolatileVolatileGasCompositeCurve::curve(const AccumulationFunction *RarAccum, const FluxFunction *RarFlux,
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
                                         int &edge){

    compositecurve.clear();
    compositecurve.type = COMPOSITE_CURVE;
    compositecurve.family = family;

    new_rarcurve.clear();
    new_rarcurve.type = RAREFACTION_CURVE;
    new_rarcurve.family = family;

    rarflux     = RarFlux;
    raraccum    = RarAccum;
    rarboundary = RarBoundary;

    // Store the first point.
    //
    add_point_to_curve(composite_initial_point, last_point_in_rarefaction, rarcurve, compositecurve);
    new_rarcurve.curve.push_back(rarcurve.curve[last_point_in_rarefaction]);

    int index_of_corresponding_point_in_rarefaction = last_point_in_rarefaction;

    RealVector old_point;
    RealVector ext_p = composite_initial_point;

    std::cout << "Max. = " << rarcurve.curve.size() << std::endl;

    for (int i = rarcurve.curve.size() - 1; i >= 0 ; i--){
        std::cout << "    i = " << i << std::endl;

        // Update (at the beginnig).
        //
        old_point = ext_p;

        int infoext = evapext->extension(rarcurve.curve[i], ext_p);

        if (infoext != EXTENSION_OK){
            std::cout << "DeadVolatileVolatileGasCompositeCurve, error while computing the extension of " << rarcurve.curve[i] << std::endl;
            return COMPOSITE_ERROR;
        }

        RealVector intersection_point;
        int info_intersect = boundary->intersection(old_point, ext_p, intersection_point, edge);

        if (info_intersect == BOUNDARY_INTERSECTION_FOUND){
            compositecurve.last_point = intersection_point;
            add_point_to_curve(intersection_point, index_of_corresponding_point_in_rarefaction, rarcurve, compositecurve);
            new_rarcurve.curve.push_back(rarcurve.curve[index_of_corresponding_point_in_rarefaction]); // TODO: WRONG!!!

            reason_why = COMPOSITE_REACHED_BOUNDARY;
            compositecurve.reason_to_stop = COMPOSITE_REACHED_BOUNDARY;

            std::cout << "Composite reached the boundary" << std::endl;

            return COMPOSITE_OK;
        }

        // Add the point.
        //
        std::cout << "Before adding the point." << std::endl;

        add_point_to_curve(ext_p, index_of_corresponding_point_in_rarefaction, rarcurve, compositecurve);
        std::cout << "After adding the point." << std::endl;

        std::cout << "Before adding the point to the rarefaction curve, index_of_corresponding_point_in_rarefaction = " << index_of_corresponding_point_in_rarefaction << std::endl;
        new_rarcurve.curve.push_back(rarcurve.curve[index_of_corresponding_point_in_rarefaction]);
        std::cout << "After adding the point to the rarefaction curve." << std::endl;

        index_of_corresponding_point_in_rarefaction--;
    }

    std::cout << "After for." << std::endl;

    // The composite curve reached the beginnig of the rarefaction.
    //
    compositecurve.final_direction = final_direction = ext_p - old_point;
    normalize(final_direction);
    normalize(compositecurve.final_direction);

    compositecurve.last_point = compositecurve.curve.back();
                    
    reason_why = COMPOSITE_COMPLETED;
    compositecurve.reason_to_stop = COMPOSITE_COMPLETED;

    return COMPOSITE_OK;
}


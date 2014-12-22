#include "CoreyQuadWaveCurveFactory.h"

CoreyQuadWaveCurveFactory::CoreyQuadWaveCurveFactory(const FluxFunction *ff, const AccumulationFunction *gg, 
                                                               const Boundary *bb, const ODE_Solver *o, 
                                                               RarefactionCurve *r, ShockCurve *s, CompositeCurve *c, 
                                                               CoreyQuadSubPhysics *tpfsp): WaveCurveFactory(gg, ff, bb, o, r, s, c), coreyquadsubphysics_(tpfsp){
}

CoreyQuadWaveCurveFactory::~CoreyQuadWaveCurveFactory(){
}

int CoreyQuadWaveCurveFactory::wavecurve(int type, const RealVector &initial_point, int family, int increase, HugoniotContinuation *h, 
                                              void *linobj, double (*linear_function)(void *o, const RealVector &p),
                                              WaveCurve &hwc, 
                                              int &wavecurve_stopped_because, int &edge){

    std::cout << "CoreyQuadWaveCurveFactory::wavecurve(), type = " << type << std::endl;

    if (type == COREYQUADWAVECURVEFACTORY_GENERIC_POINT){
        WaveCurveFactory::wavecurve(WAVECURVEFACTORY_GENERIC_POINT, initial_point, family, increase, h, linobj, linear_function, hwc, wavecurve_stopped_because, edge);
    }
    // Vertex types.
    //
    else if (type == COREYQUADWAVECURVEFACTORY_O_TO_B ||
             type == COREYQUADWAVECURVEFACTORY_O_TO_W ||
             type == COREYQUADWAVECURVEFACTORY_O_TO_G ||
             type == COREYQUADWAVECURVEFACTORY_W_TO_E ||
             type == COREYQUADWAVECURVEFACTORY_W_TO_G ||
             type == COREYQUADWAVECURVEFACTORY_W_TO_O ||
             type == COREYQUADWAVECURVEFACTORY_G_TO_D ||
             type == COREYQUADWAVECURVEFACTORY_G_TO_W ||
             type == COREYQUADWAVECURVEFACTORY_G_TO_O
            ){

        if (family != 1 || increase != RAREFACTION_SPEED_SHOULD_DECREASE) return COREYQUADWAVECURVEFACTORY_INVALID_PARAMETERS;

        RealVector bgnpoint, endpoint, dirpoint;
        
        if      (type == COREYQUADWAVECURVEFACTORY_W_TO_E ||
                 type == COREYQUADWAVECURVEFACTORY_W_TO_G ||
                 type == COREYQUADWAVECURVEFACTORY_W_TO_O){

            bgnpoint = coreyquadsubphysics_->W();
            dirpoint = coreyquadsubphysics_->E();

            if      (type == COREYQUADWAVECURVEFACTORY_W_TO_E) endpoint = coreyquadsubphysics_->E();
            else if (type == COREYQUADWAVECURVEFACTORY_W_TO_G) endpoint = coreyquadsubphysics_->G();
            else if (type == COREYQUADWAVECURVEFACTORY_W_TO_O) endpoint = coreyquadsubphysics_->O();
        }
        else if (type == COREYQUADWAVECURVEFACTORY_G_TO_D ||
                 type == COREYQUADWAVECURVEFACTORY_G_TO_W ||
                 type == COREYQUADWAVECURVEFACTORY_G_TO_O){

            bgnpoint = coreyquadsubphysics_->G();
            dirpoint = coreyquadsubphysics_->D();

            if      (type == COREYQUADWAVECURVEFACTORY_G_TO_D) endpoint = coreyquadsubphysics_->D();
            else if (type == COREYQUADWAVECURVEFACTORY_G_TO_W) endpoint = coreyquadsubphysics_->W();
            else if (type == COREYQUADWAVECURVEFACTORY_G_TO_O) endpoint = coreyquadsubphysics_->O();
        }
        else if (type == COREYQUADWAVECURVEFACTORY_O_TO_B ||
                 type == COREYQUADWAVECURVEFACTORY_O_TO_G ||
                 type == COREYQUADWAVECURVEFACTORY_O_TO_W){

            bgnpoint = coreyquadsubphysics_->O();
            dirpoint = coreyquadsubphysics_->B();

            if      (type == COREYQUADWAVECURVEFACTORY_O_TO_B) endpoint = coreyquadsubphysics_->B();
            else if (type == COREYQUADWAVECURVEFACTORY_O_TO_G) endpoint = coreyquadsubphysics_->G();
            else if (type == COREYQUADWAVECURVEFACTORY_O_TO_W) endpoint = coreyquadsubphysics_->W();
        }

        RealVector direction = endpoint - bgnpoint;
        normalize(direction);
        double epsilon = 1e-3;

        RealVector p = bgnpoint + epsilon*(dirpoint - bgnpoint);

        ReferencePoint ref(p, f, g, 0);

        hwc.family          = family;
        hwc.increase        = increase;
        hwc.reference_point = ref;

        // TODO: Add the vertex to the beginning of the first curve within hwc.
        //
        Liu_half_wavecurve(ref, p, 1 /*family*/, RAREFACTION_SPEED_SHOULD_DECREASE /*increase*/, SHOCK_CURVE, direction, linobj, linear_function, hwc, wavecurve_stopped_because, edge);
    }
    else if (type == COREYQUADWAVECURVEFACTORY_GW_SIDE ||
             type == COREYQUADWAVECURVEFACTORY_GO_SIDE ||
             type == COREYQUADWAVECURVEFACTORY_WO_SIDE
            ){
        if (family != 0 || increase != RAREFACTION_SPEED_SHOULD_INCREASE) return COREYQUADWAVECURVEFACTORY_INVALID_PARAMETERS;

        RealVector p0, p1;

        if      (type == COREYQUADWAVECURVEFACTORY_GW_SIDE){
            p0 = coreyquadsubphysics_->G();
            p1 = coreyquadsubphysics_->W();
        }
        else if (type == COREYQUADWAVECURVEFACTORY_GO_SIDE){
            p0 = coreyquadsubphysics_->G();
            p1 = coreyquadsubphysics_->O();
        }
        else if (type == COREYQUADWAVECURVEFACTORY_WO_SIDE){
            p0 = coreyquadsubphysics_->W();
            p1 = coreyquadsubphysics_->O();
        }

        RealVector p = project_point_onto_line_2D(initial_point, p0, p1);

        ReferencePoint ref(p, f, g, 0);

        RealVector direction;
        double dd;

        int info_initialize = rarefactioncurve->initialize(p, 0 /*family*/, RAREFACTION_SPEED_SHOULD_INCREASE /*increase*/, direction, dd);

        hwc.family          = family;
        hwc.increase        = increase;
        hwc.reference_point = ref;

        // TODO: Add the vertex to the beginning of the first curve within hwc.
        //
        Liu_half_wavecurve(ref, p, 0 /*family*/, RAREFACTION_SPEED_SHOULD_INCREASE /*increase*/, RAREFACTION_CURVE, direction, linobj, linear_function, hwc, wavecurve_stopped_because, edge);
    }

    return WAVECURVE_OK;
}
 
void CoreyQuadWaveCurveFactory::list_of_initial_points(std::vector<int> &type, std::vector<std::string> &name) const {
    // Generic point.
    //
    WaveCurveFactory::list_of_initial_points(type, name);

    // From O vertex.
    //
    type.push_back(COREYQUADWAVECURVEFACTORY_O_TO_B);
    name.push_back(std::string("O to B"));

    type.push_back(COREYQUADWAVECURVEFACTORY_O_TO_W);
    name.push_back(std::string("O to W"));

    type.push_back(COREYQUADWAVECURVEFACTORY_O_TO_G);
    name.push_back(std::string("O to G"));

    // From W vertex.
    //
    type.push_back(COREYQUADWAVECURVEFACTORY_W_TO_E);
    name.push_back(std::string("W to E"));

    type.push_back(COREYQUADWAVECURVEFACTORY_W_TO_G);
    name.push_back(std::string("W to G"));

    type.push_back(COREYQUADWAVECURVEFACTORY_W_TO_O);
    name.push_back(std::string("W to O"));

    // From G vertex.
    //
    type.push_back(COREYQUADWAVECURVEFACTORY_G_TO_D);
    name.push_back(std::string("G to D"));

    type.push_back(COREYQUADWAVECURVEFACTORY_G_TO_W);
    name.push_back(std::string("G to W"));

    type.push_back(COREYQUADWAVECURVEFACTORY_G_TO_O);
    name.push_back(std::string("G to O"));

    // GW side.
    //
    type.push_back(COREYQUADWAVECURVEFACTORY_GW_SIDE);
    name.push_back(std::string("GW side"));

    // GO side.
    //
    type.push_back(COREYQUADWAVECURVEFACTORY_GO_SIDE);
    name.push_back(std::string("GO side"));

    // WO side.
    //
    type.push_back(COREYQUADWAVECURVEFACTORY_WO_SIDE);
    name.push_back(std::string("WO side"));

    return;
}

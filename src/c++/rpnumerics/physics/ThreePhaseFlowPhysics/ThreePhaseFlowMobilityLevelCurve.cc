#include "ThreePhaseFlowMobilityLevelCurve.h"

void ThreePhaseFlowMobilityLevelCurve::init(){
    gv = subphysics_->gridvalues();
    mobility_on_grid.resize(gv->grid.rows(), gv->grid.cols());

    for (int i = 0; i < gv->grid.rows(); i++){
        for (int j = 0; j < gv->grid.cols(); j++){
//            if (!gv->point_inside(i, j)) continue;

            mobility_on_grid(i, j).resize(3);

            mobility_on_grid(i, j)(0) = water_mobility(this, gv->grid(i, j));
            mobility_on_grid(i, j)(1) = oil_mobility(this, gv->grid(i, j));
            mobility_on_grid(i, j)(2) = gas_mobility(this, gv->grid(i, j));

        }
    }

    return;
}

ThreePhaseFlowMobilityLevelCurve::ThreePhaseFlowMobilityLevelCurve(ThreePhaseFlowSubPhysics *s): ImplicitFunction(), subphysics_(s), mobility_(s->mobility()) {
    init();
}

ThreePhaseFlowMobilityLevelCurve::~ThreePhaseFlowMobilityLevelCurve(){
}

double ThreePhaseFlowMobilityLevelCurve::water_mobility(ThreePhaseFlowMobilityLevelCurve *obj, const RealVector &p){
    return obj->mobility_->water_mobility(p);
}

double ThreePhaseFlowMobilityLevelCurve::oil_mobility(ThreePhaseFlowMobilityLevelCurve *obj, const RealVector &p){
    return obj->mobility_->oil_mobility(p);
}

double ThreePhaseFlowMobilityLevelCurve::gas_mobility(ThreePhaseFlowMobilityLevelCurve *obj, const RealVector &p){
    return obj->mobility_->gas_mobility(p);
}

int ThreePhaseFlowMobilityLevelCurve::function_on_square(double *foncub, int i, int j){
    int is_square = gv->cell_type(i, j);
    double f_aux[4];

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
//            f_aux[l * 2 + k] = (*mobilityfunction)(this, gv->grid(i + l, j + k)) - level_;
            f_aux[l * 2 + k] = mobility_on_grid(i + l, j + k)(component_) - level_;
        }
    }

    foncub[1] = f_aux[0]; // Was: foncub[0][1]
    foncub[0] = f_aux[2]; // Was: foncub[0][0]
    foncub[3] = f_aux[1]; // Was: foncub[0][2]

    // Only useful if the cell is a square.
    //
    if (is_square == CELL_IS_SQUARE) foncub[2] = f_aux[3]; // Was: foncub[0][2]

    return 1;
}

void ThreePhaseFlowMobilityLevelCurve::curve(const RealVector &ref, int type, std::vector<RealVector> &c){
    c.clear();

    if (gv != 0){
        if      (type == WATER_MOBILITY_CURVE) mobilityfunction = &water_mobility;
        else if (type == OIL_MOBILITY_CURVE)   mobilityfunction = &oil_mobility;
        else if (type == GAS_MOBILITY_CURVE)   mobilityfunction = &gas_mobility;

        if      (type == WATER_MOBILITY_CURVE) component_ = 0;
        else if (type == OIL_MOBILITY_CURVE)   component_ = 1;
        else if (type == GAS_MOBILITY_CURVE)   component_ = 2;

        level_ = (*mobilityfunction)(this, ref);

        std::vector<RealVector> curve;
        std::vector< std::deque <RealVector> > deque_curves;
        std::vector <bool> is_circular;

        int method = SEGMENTATION_METHOD;
//        int info = ContourMethod::contour2d(this, curve, deque_curves, is_circular, method);

        double rect[4];
        rect[0] = subphysics_->boundary()->minimums()(0);
        rect[1] = subphysics_->boundary()->maximums()(0);
        rect[2] = subphysics_->boundary()->minimums()(1);
        rect[3] = subphysics_->boundary()->maximums()(1);

        int res[2] = {128, 128};

        int info = ContourMethodPure::contour2d(this, (Boundary*)subphysics_->boundary(), rect, res, curve);

        c = curve;
    }

    return;
}

void ThreePhaseFlowMobilityLevelCurve::curve(double level, int type, std::vector<RealVector> &c){
    c.clear();
    init();

    if (gv != 0){
        if      (type == WATER_MOBILITY_CURVE) mobilityfunction = &water_mobility;
        else if (type == OIL_MOBILITY_CURVE)   mobilityfunction = &oil_mobility;
        else if (type == GAS_MOBILITY_CURVE)   mobilityfunction = &gas_mobility;

        if      (type == WATER_MOBILITY_CURVE) component_ = 0;
        else if (type == OIL_MOBILITY_CURVE)   component_ = 1;
        else if (type == GAS_MOBILITY_CURVE)   component_ = 2;

        level_ = level;

        std::vector<RealVector> curve;
        std::vector< std::deque <RealVector> > deque_curves;
        std::vector <bool> is_circular;

        int method = SEGMENTATION_METHOD;
//        int info = ContourMethod::contour2d(this, curve, deque_curves, is_circular, method);

        double rect[4];
        rect[0] = subphysics_->boundary()->minimums()(0);
        rect[1] = subphysics_->boundary()->maximums()(0);
        rect[2] = subphysics_->boundary()->minimums()(1);
        rect[3] = subphysics_->boundary()->maximums()(1);

        int res[2] = {128, 128};

        int info = ContourMethodPure::contour2d(this, (Boundary*)subphysics_->boundary(), rect, res, curve);

        c = curve;
    }

    return;
}

double ThreePhaseFlowMobilityLevelCurve::level(const RealVector &ref, int type){
    if      (type == WATER_MOBILITY_CURVE) return water_mobility(this, ref);
    else if (type == OIL_MOBILITY_CURVE)   return oil_mobility(this, ref);
    else if (type == GAS_MOBILITY_CURVE)   return gas_mobility(this, ref);
}


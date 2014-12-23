#include "EquationFunctionLevelCurve.h"

EquationFunctionLevelCurve::EquationFunctionLevelCurve(const RpFunction *rpf, GridValues *g): ImplicitFunction(), function_(rpf){
    gv = g;
    level_function = &base_level_function;
}

EquationFunctionLevelCurve::~EquationFunctionLevelCurve(){
}

double EquationFunctionLevelCurve::level(const RealVector &p){
    return level(p, component_);
}

double EquationFunctionLevelCurve::level(const RealVector &p, int component){
    JetMatrix jm;

    function_->jet(p, jm, 0);

    return jm.function()(component);
}

double EquationFunctionLevelCurve::base_level_function(EquationFunctionLevelCurve *obj, const RealVector &p){
    return obj->level(p);
}

int EquationFunctionLevelCurve::function_on_square(double *foncub, int i, int j){
    int is_square = gv->cell_type(i, j);
    double f_aux[4];

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            f_aux[l * 2 + k] = (*level_function)(this, gv->grid(i + l, j + k)) - level_;
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

void EquationFunctionLevelCurve::curve(const RealVector &ref, int component, std::vector<RealVector> &c){
    if (gv != 0){
        component_ = component;
        level_ = level(ref, component_);

        std::vector<RealVector> curve;
        std::vector< std::deque <RealVector> > deque_curves;
        std::vector <bool> is_circular;

        int method = SEGMENTATION_METHOD;
        int info = ContourMethod::contour2d(this, curve, deque_curves, is_circular, method);

        c = curve;
    }

    return;
}


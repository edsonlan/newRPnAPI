#include "Coincidence_Contour.h"

Coincidence_Contour::Coincidence_Contour(const Coincidence *c) : ImplicitFunction() {
    coincidence = c;
}

Coincidence_Contour::~Coincidence_Contour(){
}

int Coincidence_Contour::saturation_on_square(Coincidence_Contour *obj, double *foncub, int i, int j){
    double f_aux[4];

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            f_aux[l * 2 + k] = obj->coincidence->lambda_s(obj->gv->grid(i + l, j + k)) - obj->saturation_level;
        }
    }

    foncub[1] = f_aux[0]; // Was: foncub[0][1]
    foncub[0] = f_aux[2]; // Was: foncub[0][0]
    foncub[3] = f_aux[1]; // Was: foncub[0][2]
    foncub[2] = f_aux[3]; // Was: foncub[0][2]

    return 1;
}

int Coincidence_Contour::evaporation_on_square(Coincidence_Contour *obj, double *foncub, int i, int j){
    double f_aux[4];

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            f_aux[l * 2 + k] = obj->coincidence->lambda_e(obj->gv->grid(i + l, j + k)) - obj->evaporation_level;
        }
    }

    foncub[1] = f_aux[0]; // Was: foncub[0][1]
    foncub[0] = f_aux[2]; // Was: foncub[0][0]
    foncub[3] = f_aux[1]; // Was: foncub[0][2]
    foncub[2] = f_aux[3]; // Was: foncub[0][2]

    return 1;
}

int Coincidence_Contour::coincidence_on_square(Coincidence_Contour *obj, double *foncub, int i, int j) {
    double f_aux[4];

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            f_aux[l * 2 + k] = obj->coincidence->lambda_diff(obj->gv->grid(i + l, j + k));
        }
    }

    foncub[1] = f_aux[0]; // Was: foncub[0][1]
    foncub[0] = f_aux[2]; // Was: foncub[0][0]
    foncub[3] = f_aux[1]; // Was: foncub[0][2]
    foncub[2] = f_aux[3]; // Was: foncub[0][2]
 
    return 1;
}

int Coincidence_Contour::curve(const FluxFunction *f, const AccumulationFunction *a,
                               GridValues &g, std::vector<RealVector> &coincidence_curve) {

    g.fill_functions_on_grid(f, a);
    gv = &g;

    coincidence_curve.clear();

    fos = &coincidence_on_square;

    int info = ContourMethod::contour2d(this, coincidence_curve);

    return info;
}

int Coincidence_Contour::characteristic_speed_curve(const FluxFunction *f, const AccumulationFunction *a, 
                                              GridValues &g, 
                                              const RealVector &p, int type, 
                                              std::vector<RealVector> &curve, double &lev){

    curve.clear();
    gv = &g;

    if (type == CHARACTERISTIC_SPEED_EVAPORATION){
        fos = &evaporation_on_square;

        lev = evaporation_level = coincidence->lambda_e(p);
    }
    else {
        fos = &saturation_on_square;

        lev = saturation_level = coincidence->lambda_s(p);
    }

    int info = ContourMethod::contour2d(this, curve);

    return info;
}


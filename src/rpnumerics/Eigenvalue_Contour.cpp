#include "Eigenvalue_Contour.h"

Eigenvalue_Contour::Eigenvalue_Contour(){
    gv = 0;

    levels.push_back(0.0); // At least this level can be contoured.
    family = 0;            // At least for this family.

    level = 0.0;
    family = 0;
}

// Find the minimum and maximum lambdas, as were computed on the grid.
//
void Eigenvalue_Contour::find_minmax_lambdas(int f, double &min, double &max){
    if (gv == 0) return;

    int r = gv->grid.rows();
    int c = gv->grid.cols();

    min = max = gv->e(0)[f].r;

    for (int i = 1; i < r*c; i++){
        if (gv->point_inside(i)){
            if (gv->e(i)[f].r > max) max = gv->e(i)[f].r;
            if (gv->e(i)[f].r < min) min = gv->e(i)[f].r;
        }
    }

    return;
}

// Set the levels as a vector with arbitrary values.
// The levels will be sorted from minimum to maximum.
//
void Eigenvalue_Contour::set_levels(int f, const std::vector<double> &l){
    if (gv == 0) return;

    family = f;

    levels.clear();
    levels.resize(l.size());

    for (int i = 0; i < l.size(); i++) levels[i] = l[i];

    std::sort(levels.begin(), levels.end());

    return;
}

// Set the levels from the minimum towards the maximum, with a
// uniform separation between levels.
//
void Eigenvalue_Contour::set_levels_from_delta(int f, double delta_l){
    if (gv == 0) return;

    family = f;
    levels.clear();

    double min_lambda, max_lambda;
    find_minmax_lambdas(family, min_lambda, max_lambda);
    
    double lambda = min_lambda;
    while (lambda < max_lambda){
        levels.push_back(lambda);
        lambda += delta_l;
    }
    
    return;
}

// Set the number of levels to be uniformly distributed.
//
void Eigenvalue_Contour::set_number_levels(int f, int n){
    if (gv == 0) return;

    family = f;
    levels.clear();

    n = max(n, 2); // At least two levels.

    double min_lambda, max_lambda;
    find_minmax_lambdas(family, min_lambda, max_lambda);

    double delta_l = (max_lambda - min_lambda)/(n - 1);
    levels.resize(n);
    for (int i = 0; i < n; i++) levels[i] = min_lambda + i*delta_l;

    return;
}

// Set the levels radiating from the level at the given point, 
// with the given distance between levels.
//
void Eigenvalue_Contour::set_levels_from_point(const FluxFunction *f, const AccumulationFunction *a, GridValues &g,
                                               int fam, const RealVector &p, double delta_l){
    g.fill_eigenpairs_on_grid(f, a);
    gv = &g;




    family = fam;

    levels.clear();

    double JF[4], JG[4];

    f->fill_with_jet(p.size(), ((RealVector)p).components(), 1, 0, JF, 0);
    a->fill_with_jet(p.size(), ((RealVector)p).components(), 1, 0, JG, 0);

    std::vector<eigenpair> e;
    Eigen::eig(p.size(), JF, JG, e);
    double level = e[family].r;

    double now_level = level - delta_l;

    double min_lambda, max_lambda;
    find_minmax_lambdas(family, min_lambda, max_lambda);

    // Downwards...
    while (now_level >= min_lambda){
        levels.push_back(now_level);
        now_level -= delta_l;
    }

    now_level = level;

    // Upwards
    while (now_level <= max_lambda){
        levels.push_back(now_level);
        now_level += delta_l;
    }

    std::sort(levels.begin(), levels.end()); printf("levels.size() = %d\n", levels.size());

    return;
}

// Set the level for a given family.
//
void Eigenvalue_Contour::set_level(double l, int f){
    level = l;
    family = f;

    return;
}

// Set the level at the given point.
//
void Eigenvalue_Contour::set_level_from_point(const FluxFunction *f, const AccumulationFunction *a,
                                              int fam, const RealVector &p){
    family = fam;


    cout<<"Tamanho do p: "<<p.size()<<endl;
    double JF[p.size()*p.size()], JG[p.size()*p.size()];

    f->fill_with_jet(p.size(), ((RealVector)p).components(), 1, 0, JF, 0);
    a->fill_with_jet(p.size(), ((RealVector)p).components(), 1, 0, JG, 0);

    std::vector<eigenpair> e;
    Eigen::eig(p.size(), JF, JG, e);
//    levels.push_back(e[family].r); //printf("Level = %f\n", e[family].r);
    level = e[family].r;

    return;
}

// Foncub.
//
int Eigenvalue_Contour::function_on_square(double *foncub, int i, int j){
    int is_square = gv->cell_type(i, j);

    if (!gv->eig_is_real(i + 1, j + 0)[family]) return 0;
    else foncub[0] = gv->e(i + 1, j + 0)[family].r - level;

    if (!gv->eig_is_real(i + 0, j + 0)[family]) return 0;
    else foncub[1] = gv->e(i + 0, j + 0)[family].r - level;

    if (!gv->eig_is_real(i + 0, j + 1)[family]) return 0;
    else foncub[3] = gv->e(i + 0, j + 1)[family].r - level;

    if (is_square == CELL_IS_SQUARE){
        if (!gv->eig_is_real(i + 1, j + 1)[family]) return 0;
        else foncub[2] = gv->e(i + 1, j + 1)[family].r - level;
    }

    return 1;
}

int Eigenvalue_Contour::curve(const FluxFunction *f, const AccumulationFunction *a, 
                              GridValues &g, 
                              std::vector< std::vector<RealVector> > &eigenvalues_curves,
                              std::vector<double> &out_levels){

    g.fill_eigenpairs_on_grid(f, a);
    gv = &g;

    eigenvalues_curves.clear();
    out_levels.clear();

    int info;

    for (int i = 0; i < levels.size(); i++){
        level = levels[i];
        std::vector<RealVector> curve;
        info = ContourMethod::contour2d(this, curve);

        if (curve.size() > 0){
            eigenvalues_curves.push_back(curve);
            out_levels.push_back(level);
        }
    }

    return info;
}

int Eigenvalue_Contour::curve(const FluxFunction *f, const AccumulationFunction *a, 
                              GridValues &g, 
                              std::vector<RealVector> &curve, double &l){

    g.fill_eigenpairs_on_grid(f, a);
    gv = &g;

    curve.clear();

    int info = ContourMethod::contour2d(this, curve);
    l = level;

    return info;
}


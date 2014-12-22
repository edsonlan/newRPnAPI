#include <iostream>

#include "Hugoniot_TP.h"

int Hugoniot_TP::classified_curve(const FluxFunction *f, const AccumulationFunction *a,
        GridValues &g, /* const RealVector &r,*/
        const ReferencePoint &r,
        std::vector<HugoniotPolyLine> &hugoniot_curve,
        std::vector<bool> &circular) {

    return classified_curve(f, a, g, r, hugoniot_curve);
}

int Hugoniot_TP::classified_curve(const FluxFunction *f, const AccumulationFunction *a,
        GridValues &g, /* const RealVector &r,*/
        const ReferencePoint &r,
        std::vector<HugoniotPolyLine> &hugoniot_curve) {

    // Compute the Hugoniot curve as usual
    //
    vector<RealVector> vrs;

    int info = curve(f, a, g, r, vrs);


    std::vector<RealVector> transition_list;

//    ReferencePoint ref(r, f, a, 0);
    ColorCurve colorCurve(*f, *a);
//    ColorCurve colorCurve(*f, *a);
    colorCurve.classify_segmented_curve(vrs, r, hugoniot_curve, transition_list);



    return info;
}

Hugoniot_TP::Hugoniot_TP(const FluxFunction *ff, const AccumulationFunction *aa, const Boundary *bb) : HugoniotCurve(ff, aa, bb), ImplicitFunction() {
    info_ = std::string("Hugoniot_TP");
}

Hugoniot_TP::~Hugoniot_TP() {
}

double Hugoniot_TP::complete_points(const RealVector &Uplus) {
    int n = Uplus.size();
    double F[n], G[n];

    RealVector p(Uplus);
    p.resize(3);
    p.component(2) = 1.0;

    ff->fill_with_jet(3, p.components(), 0, F, 0, 0);
    aa->fill_with_jet(3, p.components(), 0, G, 0, 0);

    double G1quad = G[0] - Gref.component(0);
    double G2quad = G[1] - Gref.component(1);
    double G3quad = G[2] - Gref.component(2);

    double X12minus = Fref.component(0) * G2quad - Fref.component(1) * G1quad;
    double X13minus = Fref.component(2) * G1quad - Fref.component(0) * G3quad;
    double X23minus = Fref.component(1) * G3quad - Fref.component(2) * G2quad;

    double X12plus = F[0] * G2quad - F[1] * G1quad;
    double X13plus = F[2] * G1quad - F[0] * G3quad;
    double X23plus = F[1] * G3quad - F[2] * G2quad;

    double den = X12plus * X12plus + X13plus * X13plus + X23plus * X23plus;

    // The shock speed can be evaluated as well:
    //
    //    double Y12 = Fref.component(0) * F[1] - Fref.component(1) * F[0];
    //    double Y13 = Fref.component(2) * F[0] - Fref.component(0) * F[2];
    //    double Y23 = Fref.component(1) * F[2] - Fref.component(2) * F[1];
    //
    // shock_speed = Uref.component(2)*(Y12*X12plus + Y13*X13plus + Y23*X23plus) / den;
    //
    // which must be calculated in the post-processing of the coloring.

    //    double darcy_speedplus = Uref.component(2)*(X12minus * X12plus + X13minus * X13plus + X23minus * X23plus) / den;

    double darcy_speedplus = (X12minus * X12plus + X13minus * X13plus + X23minus * X23plus) / den;

    return darcy_speedplus;
}

int Hugoniot_TP::function_on_square(double *foncub, int i, int j) {
    double dG1, dG2, dG3;
    double F1, F2, F3;

    double F10 = Fref.component(0);
    double F20 = Fref.component(1);
    double F30 = Fref.component(2);

    double f_aux[4];

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            dG1 = gv->G_on_grid(i + l, j + k).component(0) - Gref.component(0);
            dG2 = gv->G_on_grid(i + l, j + k).component(1) - Gref.component(1);
            dG3 = gv->G_on_grid(i + l, j + k).component(2) - Gref.component(2);

            F1 = gv->F_on_grid(i + l, j + k).component(0);
            F2 = gv->F_on_grid(i + l, j + k).component(1);
            F3 = gv->F_on_grid(i + l, j + k).component(2);

            f_aux[l * 2 + k] = dG1 * (F2 * F30 - F3 * F20)
                    - dG2 * (F1 * F30 - F3 * F10)
                    + dG3 * (F1 * F20 - F2 * F10);
        }
    }

    foncub[1] = f_aux[0]; // Was: foncub[0][1]
    foncub[0] = f_aux[2]; // Was: foncub[0][0]
    foncub[3] = f_aux[1]; // Was: foncub[0][2]
    foncub[2] = f_aux[3]; // Was: foncub[0][2]

    return 1;
}

// Default
int Hugoniot_TP::curve(const FluxFunction *f, const AccumulationFunction *a,
                       GridValues &g, 
                       const RealVector &r,
                       std::vector<RealVector> &hugoniot_curve){

    ReferencePoint ref(r, f, a, 0);

    return curve(f, a, g, ref, hugoniot_curve);
}

int Hugoniot_TP::curve(const FluxFunction *f, const AccumulationFunction *a,
                       GridValues &g, 
//                       const FluxFunction *ref_f, const AccumulationFunction *ref_a,
//                       const RealVector &r,
                       const ReferencePoint &r,
                       std::vector<RealVector> &hugoniot_curve) {

    printf("Entering curve\n");

    ff = f;
    aa = a;

    g.fill_functions_on_grid(f, a);

    printf("After fill_functions_on_grid\n");

    gv = &g;
//    int n = r.size();

//    Fref.resize(n);
//    Gref.resize(n);

//    Uref = r;

//    ref_f->fill_with_jet(n, Uref.components(), 0, Fref.components(), 0, 0);
//    ref_a->fill_with_jet(n, Uref.components(), 0, Gref.components(), 0, 0);

    Uref = r.point;
    Fref = r.F;
    Gref = r.G;

    printf("After filling ref\n");

    hugoniot_curve.clear();

    int info = ContourMethod::contour2d(this, hugoniot_curve);

    RealVector Utemp(2);

    // Here the third coordinate is assigned...
    for (int i = 0; i < hugoniot_curve.size(); i++) {
        Utemp = hugoniot_curve[i];

        hugoniot_curve[i].resize(3);
        hugoniot_curve[i].component(2) = complete_points(Utemp);

    }

    return info;
}

void Hugoniot_TP::map(const RealVector &r, double &func, RealVector &map_Jacobian) {

    int n = r.size() + 1;
    map_Jacobian.resize(2);

    double F[n], JF[n][n], G[n], JG[n][n];

    RealVector p(r);
    p.resize(3);
    p.component(2) = 1.0;

    f->fill_with_jet(3, p.components(), 1, F, &JF[0][0], 0);
    a->fill_with_jet(3, p.components(), 1, G, &JG[0][0], 0);

    double dG0 = G[0] - Gref[0];
    double dG1 = G[1] - Gref[1];
    double dG2 = G[2] - Gref[2];

    double F0_aux = F[1] * Fref[2] - F[2] * Fref[1];
    double F1_aux = F[0] * Fref[2] - F[2] * Fref[0];
    double F2_aux = F[0] * Fref[1] - F[1] * Fref[0];

    func = dG0 * F0_aux - dG1 * F1_aux + dG2 * F2_aux;

    map_Jacobian.component(0) = JG[0][0] * F0_aux + dG0 * (JF[1][0] * Fref[2] - JF[2][0] * Fref[1])
            - JG[1][0] * F1_aux - dG1 * (JF[0][0] * Fref[2] - JF[2][0] * Fref[0])
            + JG[2][0] * F2_aux + dG2 * (JF[0][0] * Fref[1] - JF[1][0] * Fref[0]);
    map_Jacobian.component(1) = JG[0][1] * F0_aux + dG0 * (JF[1][1] * Fref[2] - JF[2][1] * Fref[1])
            - JG[1][1] * F1_aux - dG1 * (JF[0][1] * Fref[2] - JF[2][1] * Fref[0])
            + JG[2][1] * F2_aux + dG2 * (JF[0][1] * Fref[1] - JF[1][1] * Fref[0]);
    return;
}

bool Hugoniot_TP::improvable(void) {
    return true;
}

void Hugoniot_TP::curve(const ReferencePoint &ref, int type, std::vector<Curve> &c){
    if (subphysics_ != 0) gv = subphysics_->gridvalues();

    if (gv != 0){
        reference_point = ref;

        // TODO: Maybe these four members will be eliminated. So far they are here because of the foncub.
        // If the foncub method changes, for example, by pre-computing the grid, these members would become redundant.
        //
        Fref = reference_point.F;
        Gref = reference_point.G;

        JFref = reference_point.JF;
        JGref = reference_point.JG;

        gv->fill_functions_on_grid(f, a);

        std::vector<RealVector> hugoniot_curve;
        std::vector< std::deque <RealVector> > curves;
        std::vector <bool> is_circular;

        int method = SEGMENTATION_METHOD;
        int info = ContourMethod::contour2d(this, hugoniot_curve, curves, is_circular, method);

        for (int i = 0; i < hugoniot_curve.size()/2; i++){
            hugoniot_curve[2*i].resize(3);
            hugoniot_curve[2*i + 1].resize(3);

            Curve temp;
            temp.curve.push_back(hugoniot_curve[2*i]);
            temp.curve.push_back(hugoniot_curve[2*i + 1]);

            c.push_back(temp);
        }
    }
	
    return;
}


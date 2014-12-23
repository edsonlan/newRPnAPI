#include "ImplicitHugoniotCurve.h"

ImplicitHugoniotCurve::ImplicitHugoniotCurve(const FluxFunction *ff, const AccumulationFunction *aa, const Boundary *bb) : HugoniotCurve(ff, aa, bb), ImplicitFunction() {
    method_ = IMPLICIT_HUGONIOT;

    info_ = std::string("Generic Implicit Hugoniot Curve");
}

ImplicitHugoniotCurve::~ImplicitHugoniotCurve(){
}

/**
 *  Inside the calculus of the RH condition, when the differences of the accumulation are small, it is
 *  better to approximate the flux differences by the accumulation. Notice that epsilon must be small,
 *  however large enough!
 *     The approximation is done in order to change [ df1 = f1(i,j) - f1(Uref) ] with the aid of the 
 *  approximation :
 *                     dF         dF       dU         dF      [ dG ]^(-1)
 *                    ----   =   ----  *  ----   =   ----  *  [----]
 *                     dG         dU       dG         dU      [ dU ]
 *  Here the inverse matrix [dG/dU]^(-1) is calculated in the standard form insde the [else] part.
 **/

int ImplicitHugoniotCurve::function_on_square(double *foncub, int i, int j) {
    int is_square = gv->cell_type(i, j);
    double dg1, dg2;
    double f_aux[4];
    double epsilon = 1.0e-5;

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            dg1 = gv->G_on_grid(i + l, j + k).component(0) - Gref.component(0);
            dg2 = gv->G_on_grid(i + l, j + k).component(1) - Gref.component(1);

            if (fabs(dg1) + fabs(dg2) >= epsilon) {
                double df1 = gv->F_on_grid(i + l, j + k).component(0) - Fref.component(0);
                double df2 = gv->F_on_grid(i + l, j + k).component(1) - Fref.component(1);

                f_aux[l * 2 + k] = dg2 * df1 - dg1*df2;
            } 
            else {
                // First-order expansion of F in terms of G.
                //

                double inv_det = 1.0 / (JGref(0) * JGref(3) - JGref(1) * JGref(2) );

                f_aux[l * 2 + k] = ((JFref(0) * JGref(3) - JFref(2) * JGref(1) + JFref(1) * JGref(2) - JFref(3) * JGref(0)) * dg1 * dg2 +
                        (JFref(1) * JGref(3) - JFref(3) * JGref(1)) * dg2 * dg2 +
                        (JFref(0) * JGref(2) - JFref(2) * JGref(0)) * dg1 * dg1) * inv_det;
            }
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

void ImplicitHugoniotCurve::curve(const ReferencePoint &ref, int type, std::vector<Curve> &c){
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
            Curve temp;
            temp.curve.push_back(hugoniot_curve[2*i]);
            temp.curve.push_back(hugoniot_curve[2*i + 1]);

            c.push_back(temp);
        }
    }

    return;
}

bool ImplicitHugoniotCurve::improvable(void) {
//    return true;
    return false;
}

int ImplicitHugoniotCurve::complete(const RealVector &p0, const RealVector &p1, const RealVector &p_init, RealVector &p_completed){
    Newton_Improvement newton_improver(this);
    int newton_info = newton_improver.newton(p0, p1, p_init, p_completed);

 
   return IMPROVABLE_OK;
}

double ImplicitHugoniotCurve::speed(const RealVector &F0, const RealVector &G0, const RealVector &F1, const RealVector &G1){
    double num = 0.0;
    double den = 0.0;

    for (int i = 0; i < F0.size(); i++){
        double dd = G0(i) - G1(i);

        num += (F0(i) - F1(i))*dd;
        den += dd*dd;
    }

    return num/den;
}

void ImplicitHugoniotCurve::equilibria(const std::vector<Curve> &curve, const ReferencePoint& ref, const RealVector &p, std::vector<RealVector> &ep){
    ep.clear();

    JetMatrix Fjet(ref.point.size()), Gjet(ref.point.size());
    f->jet(p, Fjet, 0);
    a->jet(p, Gjet, 0);

    RealVector F = Fjet.function();
    RealVector G = Gjet.function();

    double sigma = speed(ref.F, F, ref.G, G);

    std::vector<RealVector> vp(2);
    std::vector<double> vspeed(2);

    for (int i = 0; i < curve.size(); i++){
        for (int j = 0; j < curve[i].curve.size()/2; j++){
            for (int k = 0; k < 2; k++){
                vp[k] = curve[i].curve[2*j + k];
                f->jet(vp[k], Fjet, 0);
                a->jet(vp[k], Gjet, 0);

                vspeed[k] = speed(ref.F, Fjet.function(), ref.G, Gjet.function());
            }

            std::sort(vspeed.begin(), vspeed.end());

            if (sigma >= vspeed[0] && sigma < vspeed[1]){
                double alpha = (sigma - vspeed[1])/(vspeed[0] - vspeed[1]);

                ep.push_back(alpha*vp[0] + (1.0 - alpha)*vp[1]);
            }
        }
    }

    return;
}


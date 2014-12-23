#include "Hugoniot_Curve.h"


// This is the classified Hugoniot given by segments
//
int Hugoniot_Curve::classified_curve(const FluxFunction *f, const AccumulationFunction *a, 
                                     GridValues &g, const ReferencePoint &r, 
                                     std::vector<HugoniotPolyLine> &hugoniot_curve) {

    // The segments by a search algorithm are stored here, to be used in the ColorCurve
    std::vector<RealVector> vrs;
    std::vector<RealVector> testeTransitionalList;

    // This is auxiliary stuff (not filled, not used)
    std::vector< std::deque <RealVector> > curves;
    std::vector < bool > circular;
    int method = SEGMENTATION_METHOD;

    // Compute the Hugoniot curve as usual
    //
    int info = curve(f, a, g, r, vrs, curves, circular, method);

    //ReferencePoint ref(r, f, a, 0);

    ColorCurve colorCurve(*f, *a);
//    ColorCurve colorCurve(*f, *a);

    colorCurve.classify_segmented_curve(vrs, r, hugoniot_curve, testeTransitionalList);
    
    return info;
}

        

// This is the classified Hugoniot given by continuous curves
//
int Hugoniot_Curve::classified_curve(const FluxFunction *f, const AccumulationFunction *a, 
                                     GridValues &g, const ReferencePoint &r, 
                                     std::vector<HugoniotPolyLine> &hugoniot_curve,
                                     std::vector<bool> &circular) {

    // The continuous curve is stored by Contour2D
    std::vector< std::deque <RealVector> > curves;

    // This is auxiliary stuff (not filled, not used)
    std::vector<RealVector> vrs;
    int method = CONTINUATION_METHOD;

    // Compute the Hugoniot curve by continuation
    //
    int info = curve(f, a, g, r, vrs, curves, circular, method);
    int no_of_curves = curves.size();

    // There is a single list of Transitional points instead as one list for curves.
    std::vector<RealVector> testeTransitionalList;
    testeTransitionalList.clear();
    hugoniot_curve.clear();


//    ReferencePoint ref(r, f, a, 0);
    ColorCurve colorCurve(*f, *a);
//    ColorCurve colorCurve(*f, *a);

    for (int i = 0; i < no_of_curves; i++) {
        HugoniotPolyLine hugoniot;
        colorCurve.classify_continuous_curve(curves[i], r, hugoniot, testeTransitionalList);
        hugoniot_curve.push_back(hugoniot);
    }
    
    return info;
}

//int Hugoniot_Curve::curve(const FluxFunction *f, const AccumulationFunction *a,
//            GridValues &g, const RealVector &r,
//            std::vector<RealVector> &hugoniot_curve){

//    return 0;
//}


Hugoniot_Curve::~Hugoniot_Curve() {
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

int Hugoniot_Curve::function_on_square(double *foncub, int i, int j) {
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

//int Hugoniot_Curve::function_on_square(double *foncub, int i, int j) {
//    int is_square = gv->cell_type(i, j);
////    double dg1, dg2;
////    double f_aux[4];
////    double epsilon = 1.0e-5;

////    for (int l = 0; l < 2; l++) {
////        for (int k = 0; k < 2; k++) {
////            dg1 = gv->G_on_grid(i + l, j + k).component(0) - Gref.component(0);
////            dg2 = gv->G_on_grid(i + l, j + k).component(1) - Gref.component(1);

////            if (fabs(dg1) + fabs(dg2) >= epsilon) {
////                double df1 = gv->F_on_grid(i + l, j + k).component(0) - Fref.component(0);
////                double df2 = gv->F_on_grid(i + l, j + k).component(1) - Fref.component(1);

////                f_aux[l * 2 + k] = dg2 * df1 - dg1*df2;
////            } 
////            else {
////                // First-order expansion of F in terms of G.
////                //

////                double inv_det = 1.0 / (JGref(0) * JGref(3) - JGref(1) * JGref(2) );

////                f_aux[l * 2 + k] = ((JFref(0) * JGref(3) - JFref(2) * JGref(1) + JFref(1) * JGref(2) - JFref(3) * JGref(0)) * dg1 * dg2 +
////                        (JFref(1) * JGref(3) - JFref(3) * JGref(1)) * dg2 * dg2 +
////                        (JFref(0) * JGref(2) - JFref(2) * JGref(0)) * dg1 * dg1) * inv_det;
////            }
////        }
////    }

//    foncub[1] = hugoniot_on_grid(i, j); // f_aux[0]; // Was: foncub[0][1]
//    foncub[0] = hugoniot_on_grid(i + 1, j); //f_aux[2]; // Was: foncub[0][0]
//    foncub[3] = hugoniot_on_grid(i, j + 1); //f_aux[1]; // Was: foncub[0][2]

//    // Only useful if the cell is a square.
//    //
//    if (is_square == CELL_IS_SQUARE) foncub[2] = hugoniot_on_grid(i + 1, j + 1); //f_aux[3]; // Was: foncub[0][2]

//    return 1;
//}

//int Hugoniot_Curve::function_on_square(double *foncub, int i, int j) {
//    int is_square = gv->cell_type(i, j);

//    foncub[1] = row1[j];     // f_aux[0]; // Was: foncub[0][1]
//    foncub[0] = row2[j];     //f_aux[2]; // Was: foncub[0][0]
//    foncub[3] = row1[j + 1]; //f_aux[1]; // Was: foncub[0][2]

//    // Only useful if the cell is a square.
//    //
//    if (is_square == CELL_IS_SQUARE) foncub[2] = row2[j + 1]; //f_aux[3]; // Was: foncub[0][2]

//    return 1;
//}

void Hugoniot_Curve::fill_row(double *v, int row, int number_of_cols){
    for (int j = 0; j < number_of_cols; j++){
        if (!gv->point_inside(row, j)) continue;

        double delta_f1 = gv->F_on_grid(row, j).component(0) - Fref.component(0);
        double delta_f2 = gv->F_on_grid(row, j).component(1) - Fref.component(1);

        double delta_g1 = gv->G_on_grid(row, j).component(0) - Gref.component(0);
        double delta_g2 = gv->G_on_grid(row, j).component(1) - Gref.component(1);

        v[j] = delta_g2*delta_f1 - delta_g1*delta_f2;
    }

    return;
}

int Hugoniot_Curve::curve(const FluxFunction *f, const AccumulationFunction *a, 
                          GridValues &g, const ReferencePoint &r,
                          std::vector<RealVector> &hugoniot_curve,
                          std::vector< std::deque <RealVector> > &curves, std::vector <bool> &is_circular,
                          const int method) {

    ff = f;
    aa = a;

    g.fill_functions_on_grid(f, a);
    gv = &g;

//    int n = r.size();

//    Fref.resize(n);
//    Gref.resize(n);

//    JFref.resize(n, n); 
//    JGref.resize(n, n);

//    RealVector p(r);

//    ff->fill_with_jet(n, p.components(), 1, Fref.components(), JFref.data(), 0);
//    aa->fill_with_jet(n, p.components(), 1, Gref.components(), JGref.data(), 0);

    Fref = r.F;
    Gref = r.G;

    JFref = r.JF;
    JGref = r.JG;

    hugoniot_curve.clear();
    curves.clear();
    is_circular.clear();

    // Fill the working matrix
//    working_matrix.resize();

    // Notice that method splitts which curve is filled

    // Populate the auxiliary grid with the values of the Hugoniot.
    //
    int number_of_rows = g.grid.rows();
    int number_of_cols = g.grid.cols();

    row1 = new double[number_of_cols];
    row2 = new double[number_of_cols];

    int info;

    fill_row(row1, 0, number_of_cols);
//    for (int i = 1; i < number_of_rows; i++){
//        fill_row(row2, i, number_of_cols);

//        // Call curve2d.
//        //
//        info = ContourMethod::contour2d(this, i - 1, i, 0, number_of_cols, hugoniot_curve, curves, is_circular, method);

//        // Swap row1 and row2.
//        //
//        double *temp = row1;
//        row1 = row2;
//        row2 = temp;
//    }

//    hugoniot_on_grid.resize(g.grid.rows(), g.grid.cols());
//    for (int i = 0; i < g.grid.rows(); i++){
//        for (int j = 0; j < g.grid.cols(); j++){

//            if (!g.point_inside(i, j)) continue;

//            double delta_f1 = g.F_on_grid(i, j).component(0) - Fref.component(0);
//            double delta_f2 = g.F_on_grid(i, j).component(1) - Fref.component(1);

//            double delta_g1 = g.G_on_grid(i, j).component(0) - Gref.component(0);
//            double delta_g2 = g.G_on_grid(i, j).component(1) - Gref.component(1);

//            hugoniot_on_grid(i, j) = delta_g2*delta_f1 - delta_g1*delta_f2;
//        }
//    }

//    // Tinker around the reference point.
//    //
//    int ref_row, ref_col;
//    g.cell(r.point, ref_row, ref_col);

//    Matrix<RealVector> vertices(2, 2);
//    for (int i = 0; i < 2; i++){
//        for (int j = 0; j < 2; j++){
//            vertices(i, j) = g.grid(ref_row + i, ref_col + j);
//        }
//    }

//    int ipos, jpos;
//    double min_distance = std::numeric_limits<double>::max();
//    for (int i = 0; i < vertices.rows(); i++){
//        for (int j = 0; j < vertices.cols(); j++){
//            double n2 = norm2_squared(r.point - vertices(i, j));

//            if (n2 < min_distance){
//                min_distance = n2;
//                ipos = i;
//                jpos = j;
//            }
//        }
//    }
//    
//    // Modify only the vertex closest to the reference point.
//    // TODO: Verify if the reference point is near the boundary and if the point below lies inside or outside the domain.
//    //
//    double inv_det = 1.0/(JGref(0)*JGref(3) - JGref(1)*JGref(2));
//    double dg1 = g.G_on_grid(ref_row + ipos, ref_col + jpos).component(0) - Gref.component(0);
//    double dg2 = g.G_on_grid(ref_row + ipos, ref_col + jpos).component(1) - Gref.component(1);

//    hugoniot_on_grid(ref_row + ipos, ref_col + jpos) = ((JFref(0) * JGref(3) - JFref(2) * JGref(1) + JFref(1) * JGref(2) - JFref(3) * JGref(0)) * dg1 * dg2 +
//                                                        (JFref(1) * JGref(3) - JFref(3) * JGref(1)) * dg2 * dg2 +
//                                                        (JFref(0) * JGref(2) - JFref(2) * JGref(0)) * dg1 * dg1)*inv_det;


    info = ContourMethod::contour2d(this, hugoniot_curve, curves, is_circular, method);

    // Clean up.
    //
    delete [] row2;
    delete [] row1;

    return info;
}

void Hugoniot_Curve::map(const RealVector &r, double &f, RealVector &map_Jacobian) {
    int n = r.size();
    map_Jacobian.resize(n);

    double F[n], JF[n][n], G[n], JG[n][n];

    RealVector p(r);

    ff->fill_with_jet(n, p.components(), 1, F, &JF[0][0], 0);
    aa->fill_with_jet(n, p.components(), 1, G, &JG[0][0], 0);

    double dg1 = G[0] - Gref[0];
    double dg2 = G[1] - Gref[1];

    // TODO: The following commented lines can be uncommented in case of the Newton_Improvement do
    //       a division by zero!
    //    double epsilon = 1.0e-3;
    //    if (fabs(dg1) + fabs(dg2) >= epsilon){
    double df1 = F[0] - Fref[0];
    double df2 = F[1] - Fref[1];

    f = dg2*df1 - dg1*df2;

    map_Jacobian.component(0) = JF[0][0] * dg2 + JG[1][0] * df1 - JF[1][0] * dg1 - JG[0][0] * df2;
    map_Jacobian.component(1) = JF[0][1] * dg2 + JG[1][1] * df1 - JF[1][1] * dg1 - JG[0][1] * df2;
    //    }
    //    else {
    //        // First-order expansion of F in terms of G.
    //        //

    //        double det1 = (JFref[0]*JGref[3] - JFref[2]*JGref[1] + JFref[1]*JGref[2] - JFref[3]*JGref[0]);
    //        double det2 = (JFref[1]*JGref[3] - JFref[3]*JGref[1]);
    //        double det3 = (JFref[0]*JGref[2] - JFref[2]*JGref[0]);

    //        f = det1*dg1*dg2 + det2*dg2*dg2 + det3*dg1*dg1;

    //        map_Jacobian.component(0) = det1*(JG[0][0]*dg2 + JG[1][0]*dg1) + 2*det2*JG[0][0]*dg1 + 2*det3*JG[1][0]*dg2;
    //        map_Jacobian.component(1) = det1*(JG[0][1]*dg2 + JG[1][1]*dg1) + 2*det2*JG[0][1]*dg1 + 2*det3*JG[1][1]*dg2;
    //    }

    return;
}

int Hugoniot_Curve::complete(const RealVector &p0, const RealVector &p1, const RealVector &p_init, RealVector &p_completed){
    Newton_Improvement newton_improver(this);
    int newton_info = newton_improver.newton(p0, p1, p_init, p_completed);

    return IMPROVABLE_OK;
}

bool Hugoniot_Curve::improvable(void) {
//    return true;
    return false;
}


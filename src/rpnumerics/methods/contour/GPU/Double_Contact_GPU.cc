#include "Double_Contact_GPU.h"
#include "Contour2x2_Method.h" // Eliminate later (maybe) 

std::ostream& operator<<(std::ostream& stream, const std::vector<float> &v){
    for (int i = 0; i < v.size(); i++) stream << v[i] << " ";

    return stream;
}

bool Double_Contact_GPU::function_on_cell(double *val, int ir, int jr, int kl, int kr){
    int domain_i, domain_j;

    if ( !(gv_left->cell_is_real(ir, jr)) ) return false;

    if      (kr == 0) {domain_i = ir;     domain_j = jr;}
    else if (kr == 1) {domain_i = ir + 1; domain_j = jr;}
    else if (kr == 2) {domain_i = ir + 1; domain_j = jr + 1;}
    else if (kr == 3) {domain_i = ir;     domain_j = jr + 1;}

//    if (!gv_right->eig_is_real(domain_i, domain_j)[right_family]) return false;

    double lr  = gv_right->e(domain_i, domain_j)[right_family].r;
    double fr  = gv_right->F_on_grid(domain_i, domain_j).component(0);

    double hur = gv_right->G_on_grid(domain_i, domain_j).component(0);
    double gr  = gv_right->F_on_grid(domain_i, domain_j).component(1);
    double hvr = gv_right->G_on_grid(domain_i, domain_j).component(1);

    // Output        
    val[0] = lambda_left[kl] - lr;
    val[1] = lambda_left[kl]*(accum_left(0, kl) - hur) - (flux_left(0, kl) - fr);
    val[2] = lr*(accum_left(1, kl) - hvr) - (flux_left(1, kl) - gr);
    
    return true;
}

// Originally this function was: preplftc, and is defined in: locimp.F.
//
// ALWAYS: flux and accum are 2x4 matrices.
//         lambda is a vector with 4 elements.
//
//     3     2
//     o-----o
//     |     |
//     |     |
//     o-----o
//     0     1 
//
//     0 = (i, j),
//     1 = (i + 1, j),
//     2 = (i + 1, j + 1),
//     3 = (i, j + 1).
//
bool Double_Contact_GPU::prepare_cell(int i, int j) {
    int domain_i, domain_j;

    for (int kr = 0; kr < 4; kr++){
        if      (kr == 0) {domain_i = i;     domain_j = j;}
        else if (kr == 1) {domain_i = i + 1; domain_j = j;}
        else if (kr == 2) {domain_i = i + 1; domain_j = j + 1;}
        else if (kr == 3) {domain_i = i;     domain_j = j + 1;}

        if ( !gv_left->eig_is_real(domain_i, domain_j)[left_family] ) return false;
        lambda_left[kr]   = gv_left->e(domain_i, domain_j)[left_family].r;

        flux_left(0, kr)  = gv_left->F_on_grid(domain_i, domain_j).component(0);
        flux_left(1, kr)  = gv_left->F_on_grid(domain_i, domain_j).component(1);

        accum_left(0, kr) = gv_left->G_on_grid(domain_i, domain_j).component(0);
        accum_left(1, kr) = gv_left->G_on_grid(domain_i, domain_j).component(1);
    }

    return true;
}

void Double_Contact_GPU::curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues *lg, int lfam,
                           const FluxFunction *rf, const AccumulationFunction *ra, GridValues *rg, int rfam,
                           std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve){
    lff = lf;
    laa = la;
    gv_left = lg;
    left_family = lfam;

    rff = rf;
    raa = ra;
    gv_right = rg;
    right_family = rfam;

    // It is assumed that the grid_value must be the same
    singular = ( (left_family == right_family) && (lg == rg) );

    gv_left->fill_eigenpairs_on_grid(lff, laa);
    gv_right->fill_eigenpairs_on_grid(rff, raa);

    left_curve.clear(); 
    right_curve.clear(); 

    Contour2x2_Method::curve2x2(this, left_curve, right_curve);

    return;
}

void Double_Contact_GPU::gpu_curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues *lg, int lfam,
                               const FluxFunction *rf, const AccumulationFunction *ra, GridValues *rg, int rfam,
                               std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve){
    lff = lf;
    laa = la;
    gv_left = lg;
    left_family = lfam;

    rff = rf;
    raa = ra;
    gv_right = rg;
    right_family = rfam;

    // It is assumed that the grid_value must be the same
    singular = ( (left_family == right_family) && (lg == rg) );

    gv_left->fill_eigenpairs_on_grid(lff, laa);
    gv_right->fill_eigenpairs_on_grid(rff, raa);

    left_curve.clear(); 
    right_curve.clear(); 

    // TODO: If the families are equal AND the grids are also equal these blocks
    //       below could be unified.

    // Left.
    //
    int left_rows = gv_left->grid.rows();
    int left_cols = gv_left->grid.cols();
    int left_size = left_rows*left_cols;

    // TODO: These are to be declared in the header.
    //
    gpu_lambda_left.resize(left_rows, left_cols); // Matrix<float>
    gpu_F_left.resize(left_rows, left_cols);      // Matrix<std::vector<float> >
    gpu_G_left.resize(left_rows, left_cols);      // Matrix<std::vector<float> >

    // TODO: This copies from the main memory to the GPU's.
    //
    for (int i = 0; i < left_rows; i++){
    for (int j = 0; j < left_cols; j++){
//        std::cout << "i = " << i << ", j = " << j << ", cell type = " << gv_left->cell_type(i, j) << std::endl;
//        if (gv_left->cell_type(i) != CELL_IS_SQUARE) continue;

//        std::cout << "    gv_left->e(" << i << ", " << j << ").size() = " << gv_left->e(i, j).size() << std::endl;

        if (gv_left->e(i, j).size() > 0) gpu_lambda_left(i) = gv_left->e(i, j)[lfam].r;

        int n = gv_left->F_on_grid(i, j).size();
        gpu_F_left(i, j).resize(n);
        gpu_G_left(i, j).resize(n);

        for (int k = 0; k < n; k++){
            gpu_F_left(i, j)[k] = gv_left->F_on_grid(i, j)(k);
            gpu_G_left(i, j)[k] = gv_left->G_on_grid(i, j)(k);
        }
    }
    }

    // Right.
    //
    int right_rows = gv_right->e.rows();
    int right_cols = gv_right->e.cols();
    int right_size = right_rows*right_cols;

    // TODO: These are to be declared in the header.
    //
    gpu_lambda_right.resize(right_rows, right_cols);
    gpu_F_right.resize(right_rows, right_cols);
    gpu_G_right.resize(right_rows, right_cols);

    // TODO: This copies from the main memory to the GPU's.
    //
    for (int i = 0; i < right_rows; i++){
    for (int j = 0; j < right_cols; j++){
//        std::cout << "i = " << i << ", j = " << j << ", cell type = " << gv_right->cell_type(i, j) << std::endl;
//        if (gv_left->cell_type(i) != CELL_IS_SQUARE) continue;

//        std::cout << "    gv_right->e(" << i << ", " << j << ").size() = " << gv_right->e(i, j).size() << std::endl;

        if (gv_right->e(i, j).size() > 0) gpu_lambda_right(i) = gv_right->e(i, j)[rfam].r;

        int n = gv_right->F_on_grid(i, j).size();
        gpu_F_right(i, j).resize(n);
        gpu_G_right(i, j).resize(n);

        for (int k = 0; k < n; k++){
            gpu_F_right(i, j)[k] = gv_right->F_on_grid(i, j)(k);
            gpu_G_right(i, j)[k] = gv_right->G_on_grid(i, j)(k);
        }
    }
    }

    // The number of cells must be a power of 2.
    // Therefore, the number of vertices must be that same power of 2 plus 1.
    //
    int n = 1;

//    int length_left = ((left_rows)*(left_rows + 1)/2)/n;
//    std::cout << "length_left = " << length_left << std::endl;

    int length_left = (left_rows)*(left_cols);
    std::cout << "length_left = " << length_left << std::endl;

    Matrix<float> val0, val1, val2;
//    gpu_allocate_arrays_in_triangle(this, length_left, val0, val1, val2);
    gpu_allocate_arrays_in_rectangle(this, length_left, val0, val1, val2);

    // These will store the indexes of the cells where changes in sign were detected.
    //
    std::vector<short int> index_l, index_r;
    gpu_fill_arrays(0, length_left, val0, val1, val2);
    gpu_check_sign(0, length_left, val0, val1, val2, index_l, index_r);

//    std::cout << "val0 = " << std::endl << val0 << std::endl << std::endl;
//    std::cout << "val1 = " << std::endl << val1 << std::endl << std::endl;
//    std::cout << "val2 = " << std::endl << val2 << std::endl << std::endl;

    std::cout << "index_l.size() = " << index_l.size() << ", index_r.size() = " << index_r.size() << std::endl;
//    for (int i = 0; i < index_l.size(); i++){
//        std::cout << "index_l = " << index_l[i] << ", i = " << index_l[i]/gv_left->grid.cols() << ", j = " << index_l[i] - (index_l[i]/gv_left->grid.cols()) << std::endl;
//        std::cout << "index_r = " << index_r[i] << ", i = " << index_r[i]/gv_right->grid.cols() << ", j = " << index_r[i] - (index_r[i]/gv_right->grid.cols()) << std::endl;
//    }

//    for (int i = 0; i < n + 1; i++){
//        std::cout << "i = " << i << std::endl;

//        int min_left = i*length_left;
//        int max_left = (i + 1)*length_left - 1;
//        gpu_fill_arrays(min_left, max_left, val0, val1, val2);

//        std::cout << "    gpu_fill_arrays passed." << std::endl;


////        // TODO: Here comes the method that checks the signs in the cells and returns a list of the viable cells.
////        //       Afterwards, the modified contour2x2 should be called.
////        //
////        gpu_check_sign(min_left, max_left, val0, val1, val2, index_l, index_r);

////        std::cout << "    gpu_check_sign passed." << std::endl;
//    }

//    return;

    Contour2x2_Method_GPU::cpu_curve2x2(this, index_l, index_r, left_curve, right_curve);

    return;
}

void Double_Contact_GPU::gpu_allocate_arrays_in_rectangle(Double_Contact_GPU *dc, int length_left, Matrix<float> &val0, Matrix<float> &val1, Matrix<float> &val2){
    int right_rows = dc->gv_right->grid.rows();
    int right_cols = dc->gv_right->grid.cols();
    int right_size = right_rows*right_cols;

    int total_size = length_left*right_size;

    val0.resize(length_left, right_size);
    val1.resize(length_left, right_size);
    val2.resize(length_left, right_size);

    return;
}

void Double_Contact_GPU::gpu_allocate_arrays_in_triangle(Double_Contact_GPU *dc, int length_left, Matrix<float> &val0, Matrix<float> &val1, Matrix<float> &val2){
    int right_rows = dc->gv_right->grid.rows();
    int right_size = right_rows*(right_rows + 1)/2;

    int total_size = length_left*right_size;

    val0.resize(length_left, right_size);
    val1.resize(length_left, right_size);
    val2.resize(length_left, right_size);

    return;
}

// TODO: See if it is necessary to store the indices of the elements at left and right.
//       This coming and going makes it hard on the parallelization.
//
void Double_Contact_GPU::gpu_fill_arrays(int min_left, int max_left, Matrix<float> &val0, Matrix<float> &val1, Matrix<float> &val2){
    // TODO: These are computed here, but they could come from the outside.
    //
    int min_right, max_right;
    min_right = 0;
    max_right = val0.cols();  // Only works OK with triangles or squares.
//    std::cout << "val0.rows() = " << val0.rows() << ", val0.cols() = " << val0.cols() << ", max_right = " << max_right << std::endl;
//    std::cout << "gpu_lambda_left = " << gpu_lambda_left.rows() << " x " << gpu_lambda_left.cols() << std::endl;
//    std::cout << "gpu_F_left = " << gpu_F_left.rows() << " x " << gpu_F_left.cols() << std::endl;
//    std::cout << "gpu_G_left = " << gpu_G_left.rows() << " x " << gpu_G_left.cols() << std::endl;
//    std::cout << "min_left = " << min_left << ", max_left = " << max_left << std::endl;

    float nl = (float)gv_left->grid.rows() + 1e-3;
    float deltaxl = gv_left->grid_resolution(0);
    float deltayl = gv_left->grid_resolution(1);

    float deltaxr = gv_right->grid_resolution(0);

    // TODO: Initial and end values must be flexible to use parallelism.
    //
    for (int l = min_left; l < max_left; l++){
        // These variables are explicitly used here
        // only to make the code clear.
        //
        int l_cols = gpu_lambda_left.cols();
        int il = std::floor(l/l_cols);
        int jl = l - il*l_cols;

//        std::cout << "l = " << l << ", l_cols = " << l_cols << std::endl;
//        std::cout << "    il = " << il << ", jl = " << jl << std::endl;

        float lambda_left = gpu_lambda_left(il, jl);      // Lambda
//        std::cout << "    lambda_left = " << lambda_left << std::endl;

        std::vector<float> F_left = gpu_F_left(il, jl);         // Flux
//        for (int m = 0; m < F_left.size(); m++) std::cout << "    F_left[" << m << "] = " << F_left[m] << std::endl;

        std::vector<float> G_left = gpu_G_left(il, jl); // Accumulation
//        for (int m = 0; m < G_left.size(); m++) std::cout << "    G_left[" << m << "] = " << G_left[m] << std::endl;

//        std::cout << "    here" << std::endl;

        for (int r = min_right; r < max_right; r++){
            // These variables are explicitly used here
            // only to make the code clear.
            //
            int r_cols = gpu_lambda_right.cols();
            int ir = std::floor(r/r_cols);
            int jr = r - ir*r_cols;

//            std::cout << "    r = " << r << ", r_cols = " << r_cols << ", max_right = " << max_right << std::endl;
//            std::cout << "        ir = " << ir << ", jr = " << jr << std::endl;

            float lambda_right = gpu_lambda_right(ir, jr);      // Lambda
//            std::cout << "        lambda_right = " << lambda_right << std::endl;

//            std::cout << "        gpu_F_right = " << gpu_F_right.rows() << " x " << gpu_F_right.cols() << std::endl;
            std::vector<float> F_right = gpu_F_right(ir, jr);         // Flux
//            for (int m = 0; m < F_right.size(); m++) std::cout << "    F_right[" << m << "] = " << F_right[m] << std::endl;

            std::vector<float> G_right = gpu_G_right(ir, jr); // Accumulation
//            for (int m = 0; m < G_right.size(); m++) std::cout << "    G_right[" << m << "] = " << G_right[m] << std::endl;

//            val0(l, r) = lambda_left - lambda_right; // (lambda_left - lambda_right > 0.0 ? 1 : -1); Try using a Matrix<short int>.
            val0(l, r) = (nl - 2.0*il)*deltaxl; // (lambda_left - lambda_right > 0.0 ? 1 : -1); Try using a Matrix<short int>.
//            std::cout << "        =>>> val0(" << l << ", " << r << ") = " << val0(l, r) << std::endl;

//            std::cout << "lambda_left = " << lambda_left << std::endl;
//            std::cout << "G_left[0]  = " << G_left[0] << std::endl;
//            std::cout << "G_right[0] = " << G_right[0] << std::endl;
//            std::cout << "F_left[0]  = " << F_left[0] << std::endl;
//            std::cout << "F_right[0] = " << F_right[0] << std::endl;

//            val1(l, r) = lambda_left *(G_left[0] - G_right[0]) - (F_left[0] - F_right[0]);
            val1(l, r) = (nl - 2.0*jl)*deltayl;
//            std::cout << "        =>>> val1(" << l << ", " << r << ") = " << val1(l, r) << std::endl;

//            val2(l, r) = lambda_right*(G_left[1] - G_right[1]) - (F_left[1] - F_right[1]);
            val2(l, r) = (nl - 2.0*ir)*deltaxr;

//            std::cout << "l = " << l << ", r = " << r << ", ir = " << ir << ", deltaxr = " << deltaxr << ", (nl - 2.0*ir)*deltaxr = " << (nl - 2.0*ir)*deltaxr << std::endl;
//            std::cout << "        =>>> val2(" << l << ", " << r << ") = " << val2(l, r) << std::endl;
        }
    }

    std::cout << "Done with the array-filling part." << std::endl;

    return;
}

void Double_Contact_GPU::gpu_check_sign(int min_left, int max_left, const Matrix<float> &val0, const Matrix<float> &val1, const Matrix<float> &val2, 
                                        std::vector<short int> &index_l, std::vector<short int> &index_r){
    int min_right, max_right;
    min_right = 0;
    max_right = val0.cols();

    // TODO: Initial and end values must be flexible to use parallelism.
    //

    // The values of the functions are needed in the four vertices of every cell.
    // The cells are being indexed by a single index:
    //
    //     I = i*cols + j,
    //
    // where i, j are the indices of the lower-left vertex of the abovementioned cell.
    //
    // The indices of the rest of the vertices, that form the cell are, in the usual two-index fashion:
    //
    //     (i, j+1), (i+1, j), (i+1, j+1).
    //
    // These need to be translated to the single-index fashion, thus:
    //
    // (i, j+1)   ==> I + 1,
    // (i+1, j)   ==> I + cols and
    // (i+1, j+1) ==> I + cols + 1.
    //
    // This applies to both indices of val0, etc.
    //

    int Lcols = gv_left->grid.cols();
    int Rcols = gv_right->grid.cols();

    for (int l = min_left; l < max_left; l++){
//        int l0, l1, l2, l3; // Counterclockwise.
//        l0 = l;
//        l1 = l + 1;
//        l3 = l + Lcols;
//        l2 = l3 + 1;

        std::vector<unsigned short int> lindex(4);
        lindex[0] = l;
        lindex[1] = l + 1;
        lindex[3] = l + Lcols;
        lindex[2] = lindex[3] + 1;

        for (int r = min_right; r < max_right; r++){
//            unsigned short int r0, r1, r2, r3; // Counterclockwise.
//            r0 = r;
//            r1 = r + 1;
//            r3 = r + Rcols;
//            r2 = r3 + 1;

            std::vector<unsigned short int> rindex(4);
            rindex[0] = r;
            rindex[1] = r + 1;
            rindex[3] = r + Rcols;
            rindex[2] = rindex[3] + 1;

            bool val0_change = false, val1_change = false, val2_change = false;

            // Check the first element against the rest of the elements of the first column...
            //
            float ref_val0 = val0(l, r);
            for (int i = 1; i < 4; i++){
                if (ref_val0*val0(i, r) < 0.0){
                    val0_change = true;
                    break;
                }
            }

            // ...and if nothing was found try against all the other rows & columns.
            //
            if (!val0_change){
                for (int i = 0; i < 4; i++){
                for (int j = 1; j < 4; j++){
                    if (ref_val0*val0(lindex[i], rindex[j]) < 0.0){
                        val0_change = true;
                        break;
                    }
                }
                }

                if (val0_change){
                    // Check//            // These variables are explicitly used here
//            // only to make the code clear.
//            //
//            int r_cols = gpu_lambda_left.cols();
//            int ir = std::floor(r/r_cols);
//            int jr = r - ir*r_cols;

//            std::cout << "    r = " << r << ", r_cols = " << r_cols << std::endl;
//            std::cout << "        ir = " << ir << ", jr = " << jr << std::endl;

//            float lambda_right = gpu_lambda_right(ir, jr);      // Lambda
//            std::cout << "        lambda_right = " << lambda_right << std::endl;

//            std::cout << "        gpu_F_right = " << gpu_F_right.rows() << " x " << gpu_F_right.cols() << std::endl;
//            std::vector<float> F_right = gpu_F_right(ir, jr);         // Flux
//            for (int m = 0; m < F_right.size(); m++) std::cout << "    F_right[" << m << "] = " << F_right[m] << std::endl;

//            std::vector<float> G_right = gpu_G_right(ir, jr); // Accumulation
//            for (int m = 0; m < G_right.size(); m++) std::cout << "    G_right[" << m << "] = " << G_right[m] << std::endl; the first element against the rest of the elements of the first column...
                    //
                    float ref_val1 = val1(l, r);
                    for (int i = 1; i < 4; i++){
                        if (ref_val1*val1(i, r) < 0.0){
                            val1_change = true;
                            break;
                        }
                    }

                    // ...and if nothing was found try against all the other rows & columns.
                    //
                    if (!val1_change){
                        for (int i = 0; i < 4; i++){
                        for (int j = 1; j < 4; j++){
                            if (ref_val1*val1(lindex[i], rindex[j]) < 0.0){
                                val1_change = true;
                                break;
                            }
                        }
                        }
                    }

                    if (val1_change){
                        float ref_val2 = val2(l, r);
                        for (int i = 1; i < 4; i++){
                            if (ref_val2*val2(i, r) < 0.0){
                                val2_change = true;
                                break;
                            }
                        }

                        // ...and if nothing was found try against all the other rows & columns.
                        //
                        if (!val2_change){
                            for (int i = 0; i < 4; i++){
                            for (int j = 1; j < 4; j++){
                                if (ref_val2*val2(lindex[i], rindex[j]) < 0.0){
                                    val2_change = true;
                                    break;
                                }
                            }
                            }

                            // If the sign of all of the three functions changed, add this cell to
                            // the list where Contour2x2 will work.
                            //
                            if (val2_change){
                                index_l.push_back(l);
                                index_r.push_back(r);
                            }

                        }
                    } // For val2
                } // For val1
            } // For val0



        }
    }

    return;
}

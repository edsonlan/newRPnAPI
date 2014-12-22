#include "eigen.h"

#include "FluxFunction.h"
#include "AccumulationFunction.h"

// Initialize the value of the epsilon
double Eigen::epsilon(1e-10);

// Function to compare the ratio alphar/beta, to be passed to the std::algorithm::sort method.
// TODO: Check this one.
bool Eigen::eigen_compare(const eigenpair &n1, const eigenpair &n2){
    return n1.r < n2.r;
}

// Function to compare the ratio alphar/beta, to be passed to the std::algorithm::sort method.
// TODO: Check this one.

//// TPCW!!!
//// In this case the eigenpairs are not sorted according to their eigenvalues,
//// but according to their eigenvectors (if they are horizontal or not).
////
//bool Eigen::eigen_compare(const eigenpair &n1, const eigenpair &n2){
//    RealVector r1(2), r2(2);
//    for (int i = 0; i < 2; i++){
//        r1(i) = n1.vrr[i];
//        r2(i) = n2.vrr[i];
//    }

//    normalize(r1);
//    normalize(r2);

//    return std::abs(r1(1)) < std::abs(r2(1));
//}

bool Eigen::eigen_compare_using_eigenvalues(const eigenpair &n1, const eigenpair &n2){
    return n1.r < n2.r;
}

bool Eigen::eigen_compare_using_eigenvectors(const eigenpair &n1, const eigenpair &n2){
    RealVector r1(2), r2(2);
    for (int i = 0; i < 2; i++){
        r1(i) = n1.vrr[i];
        r2(i) = n2.vrr[i];
    }

    normalize(r1);
    normalize(r2);

    return std::abs(r1(1)) < std::abs(r2(1));
}

// Function to be used to compare the eigenpairs.
//
bool (*Eigen::eigen_sort_function)(const eigenpair&, const eigenpair&) = &eigen_compare_using_eigenvalues;

// Transpose a square matrix, rewriting the original.
//
void Eigen::transpose(int n, double *A){
    int i, j;
    double temp;
    
    for (i = 0; i < n; i++){
        for (j = i; j < n; j++){
            temp = A[i*n + j];
            A[i*n + j] = A[j*n + i];
            A[j*n + i] = temp;
        }
    }
    return;
}

// To fill an array of eigen structs. It is assumed that the matrix containing the
// eigenvectors returned by Lapack was tranposed before being fed to
// this function.
void Eigen::fill_eigen(int n, struct eigenpair e[], double rp[], double ip[], double vl[], double vr[]){
    int i, j;
     
    for (i = 0; i < n; i++){
        if (fabs(ip[i]) < Eigen::epsilon){ // Eigenvalue is real
            e[i].r = rp[i];
            e[i].i = 0;
            for (j = 0; j < n; j++){
                e[i].vlr[j] = vl[j*n + i];
                e[i].vli[j] = 0;
                e[i].vrr[j] = vr[j*n + i];
                e[i].vri[j] = 0;                
            }
        }
        else{                   // Eigenvalue is complex
            if (ip[i] > 0){     // Make sure the eigenvalue's imaginary part is positive.
                                // In this case the i-th column of v contains the real part
                                // of the eigenvector and the (i + 1)-th column contains the
                                // imaginary part of the eigenvector. 
                               
                // If the eigenvalues are complex they are returned by Lapack in conjugated pairs,
                // the first of the pair having positive imaginary part
                e[i].r = rp[i];
                e[i + 1].r = rp[i];
/*                
                // Begin of [1] 
                // This code replaces the block below marked as [2].
                // See also blocks [3] and [4].
                // It should be better than the former, since it thresholds very small values of the
                // imaginary part of the eigenvector.
                double temp_i = ip[i];
                if (fabs(temp_i) < EPS) temp_i = 0;
                e[i].i = temp_i;
                e[i + 1].i = -temp_i;
                // End of [1]
*/                
                
                // Begin of [2]
                // This code was replaced by the code above, marked as [1].
                // See also blocks [3] and [4].
                e[i].i = ip[i];
                e[i + 1].i = -ip[i];
                // End of [2] 
                
                
/*
                // Begin of [3] 
                // This code replaces the block below marked as [4].
                // See also blocks [1] and [2].
                // It should be better than the former, since it thresholds very small values of the
                // imaginary part of the eigenvector.
                double temp_vi;
                for (j = 0; j < n; j++){
                    // Real part of the pair of complex eigenvectors:
                    e[i].vr[j]     = v[j*n + i];
                    e[i + 1].vr[j] = v[j*n + i];
                    
                    // For the imaginary part of the pair of complex eigenvectors:
                    temp_vi = v[j*n + i + 1];
                    if (fabs(temp_vi) < EPS) temp_vi = 0;
                    e[i].vi[j]     = temp_vi;
                    e[i + 1].vi[j] = -temp_vi;
                }
*/
                
                // Begin of [4]
                // This code was replaced by the code above, marked as [3].
                // See also blocks [1] and [2].
                for (j = 0; j < n; j++){
                    e[i].vlr[j] = vl[j*n + i];
                    e[i].vli[j] = vl[j*n + i + 1];
                    
                    e[i + 1].vlr[j] = vl[j*n + i];
                    e[i + 1].vli[j] = -vl[j*n + i + 1];
                    
                    e[i].vrr[j] = vr[j*n + i];
                    e[i].vri[j] = vr[j*n + i + 1];
                    
                    e[i + 1].vrr[j] = vr[j*n + i];
                    e[i + 1].vri[j] = -vr[j*n + i + 1];                    
                }
                // End of [4]
                
                
                i++;
            }
            else{              // This should never happen, but just in case...
                printf("Problem in fill_eigen! i = %d\n", i);
            }
        }
    }    
    return;
}

// Engine of sort_eigen (below)
int Eigen::eigen_comp(const void *p1, const void *p2){
   eigenpair *sp1 = (eigenpair*)p1;
   eigenpair *sp2 = (eigenpair*)p2;

   double temp = sp1->r - sp2->r;

   if (temp > 0) return 1;
   else if (temp < 0) return -1;
   else return 0;
}

// Sort e[] according to the real part of the eigenvalues
void Eigen::sort_eigen(int n, eigenpair *e){
    qsort(e, n, sizeof(e[0]), eigen_comp);
    return;
}

// Compute all the eigenpairs of a matrix and fill a
// vector of said structs.
//
//      n: Dimension of the space
//      A: Square matrix
//     ve: Vector where the eigenpairs will be stored.
//
// This function uses LAPACK's dgeev. Eigenvectors will be
// normalized.
//
int Eigen::eig(int n, const double *A, vector<eigenpair> &ve){
    int lda = n, lwork = 5*n, ldvr = n, ldvl = n;
    int i, j, info;
    double vr[n][n], vl[n][n];
    double work[5*n], wi[n], wr[n];

    // Create a transposed copy of A to be used by LAPACK's dgeev:
    double B[n][n];
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++) B[j][i] = A[i*n + j];
    }

    dgeev_("V", "V", &n, &B[0][0], &lda, &wr[0], &wi[0], 
           &vl[0][0], &ldvl, &vr[0][0], &ldvr, &work[0], &lwork, 
           &info);

    // Process the results

    if (info != 0) return info;
    else {
        transpose(n, &vl[0][0]); // ...or else...
        transpose(n, &vr[0][0]); // ...or else...

        eigenpair e[n];
        for (int i = 0; i < n; i++){
            e[i].vlr.resize(n);
            e[i].vli.resize(n);
            e[i].vrr.resize(n);
            e[i].vri.resize(n);
        }

        fill_eigen(n, e, &wr[0], &wi[0], &vl[0][0], &vr[0][0]);
        sort_eigen(n, e);

        ve.clear();
        ve.resize(n);
        for (int i = 0; i < n; i++) ve[i] = e[i];
        return 0;// SUCCESSFUL_PROCEDURE;

    }

    return info;
}

// Generalized eigenproblem
// TODO: Add an option to compute only the eigenvalues (no eigenvectors).
int Eigen::eig(int n, const double *A, const double *B, vector<eigenpair> &vge){
    // Dimension
    int dim = n;

    // Copy A and B (transposed)
    double AA[dim*dim], BB[dim*dim];
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            AA[i*dim + j] = A[j*dim + i];
            BB[i*dim + j] = B[j*dim + i];
        }
    }

    // Some parameters
    int lda = max(1, n);
    int ldb = max(1, n);

    // Eigenvalues of the form:
    //
    // (alphar[j] + i*alphai[j])/beta[j] for j = 0,..., (n - 1).
    //
    double alphar[dim], alphai[dim], beta[dim];

    // The leading dimension of the matrix vl. 
    // ldlvl >= 1, and if jovl = "V", ldvl >= n.
    // 
    int ldvl = n;

    double vl[ldvl*dim]; 
    
    // The leading dimension of the matrix vr. 
    // ldvr >= 1, and if jobvr = "V", ldvr >= n.
    int ldvr = n;

    double vr[ldvr*dim]; 

    // Working array
    int lwork = max(1, 8*dim);
    double work[lwork];

    // Info
    int info = 0;

    // Invoke LAPACK
    dggev_("V", "V", &dim, 
           AA, &lda, 
           BB, &ldb, 
           alphar, alphai, beta, 
           vl, &ldvl, 
           vr, &ldvr,
           work, &lwork,
           &info);

    transpose(dim, vl);
    transpose(dim, vr);

    // Success!
    if (info == 0) {
//        for (int i = 0; i < dim; i++) printf("alphar[%d] = %g, alphai[%d] = %g, beta[%d] = %g\n", i, alphar[i], i, alphai[i], i, beta[i]);

        // Abort if some beta is smaller than sum_i(abs(alphar(i)) + abs(alphai(i))).
        double sum = 0;
        for (int i = 0; i < dim; i++) sum += fabs(alphar[i]) + fabs(alphai[i]);
        sum *= epsilon;

        // Quick and dirty
        int max = 0;
        for (int i = 0; i < dim; i++){
            if (fabs(beta[i]) > sum) max++; //return _EIGEN_BETA_NEAR_ZERO_;
        }
        if (max == 0) return _EIGEN_BETA_NEAR_ZERO_;

        // Otherwise, carry on
        vge.clear();
        vge.resize(max);

        for (int i = 0; i < max; i++){
            vge[i].vlr.resize(dim);
            vge[i].vli.resize(dim);
            vge[i].vrr.resize(dim);
            vge[i].vri.resize(dim);
        }

        int pos = 0;
        for (int i = 0; i < dim; i++){
            if (fabs(beta[i]) > sum){
                if (fabs(alphai[i]) < epsilon){
                    // Eigenvalue is real.

                    vge[pos].r = alphar[i]/beta[i];
                    vge[pos].i = alphai[i]/beta[i];

                    // Normalize
                    double norm_l(0), norm_r(0);
                    for (int j = 0; j < dim; j++){
                        norm_l += vl[i + j*dim]*vl[i + j*dim];
                        norm_r += vr[i + j*dim]*vr[i + j*dim];
                    }
                    norm_l = sqrt(norm_l);
                    norm_r = sqrt(norm_r);

                    for (int j = 0; j < dim; j++){
                        vge[pos].vlr[j] = vl[i + j*dim]/norm_l;
                        vge[pos].vli[j] = 0;

                        vge[pos].vrr[j] = vr[i + j*dim]/norm_r;
                        vge[pos].vri[j] = 0;
                    }
                }
                else {
                    // Eigenvalue is complex.
                    // In this case the i-th column of v contains the real part
                    // of the eigenvector and the (i + 1)-th column contains the
                    // imaginary part of the eigenvector. 
                                   
                    // If the eigenvalues are complex they are returned by Lapack in conjugated pairs,
                    // the first of the pair having positive imaginary part.

                    vge[pos].r = alphar[i]/beta[i];
                    vge[pos].i = alphai[i]/beta[i];
    
                    vge[pos + 1].r = alphar[i]/beta[i];
                    vge[pos + 1].i = -alphai[i]/beta[i];

                    // Normalize
                    double norm_l(0), norm_r(0);
                    for (int j = 0; j < dim; j++){
                        norm_l += vl[i + j*dim]*vl[i + j*dim] + vl[i + 1 + j*dim]*vl[i + 1 + j*dim];
                        norm_r += vr[i + j*dim]*vr[i + j*dim] + vr[i + 1 + j*dim]*vr[i + 1 + j*dim];
                    }
                    norm_l = sqrt(norm_l);
                    norm_r = sqrt(norm_r);

                    for (int j = 0; j < dim; j++){
                        vge[pos].vlr[j]     = vl[i + j*dim]/norm_l;
                        vge[pos].vli[j]     = vl[i + 1 + j*dim]/norm_l;
    
                        vge[pos + 1].vlr[j] = vl[i + j*dim]/norm_l;
                        vge[pos + 1].vli[j] = -vl[i + 1 + j*dim]/norm_l;

                        vge[pos].vrr[j]     = vr[i + j*dim]/norm_r;
                        vge[pos].vri[j]     = vr[i + 1 + j*dim]/norm_r;

                        vge[pos + 1].vrr[j] = vr[i + j*dim]/norm_r;
                        vge[pos + 1].vri[j] = -vr[i + 1 + j*dim]/norm_r;
                    }

                    // Skip one step!
                    i++;
                }
                pos++;
            }
            else {
                //printf("Eigenvalue discarded: %d\n", pos);
            }
        }

//        sort(vge.begin(), vge.end(), eigen_compare);
        sort(vge.begin(), vge.end(), eigen_sort_function);

    }

    return info;
}

// Generalized eigenproblem for a given family
int Eigen::eig(int n, const double *A, const double *B, int family, eigenpair &ep){
    std::vector<eigenpair> e;
    int info = eig(n, A, B, e);
    
//    std::cout << "Eigen, new method. info = " << info << std::endl;
    print_eigen(e);

    if (info == 0) ep = e[family];
    
    return info;
}

// Generalized eigenproblem for a given family, returning ONLY the right-eigenvector
int Eigen::eig(int n, const double *A, const double *B, int family, RealVector &r){
    eigenpair e;
    int info = eig(n, A, B, family, e);
    
    std::vector<eigenpair> ee;
    ee.push_back(e);
    print_eigen(ee);
    
    if (info == 0){
        r.resize(n);
        for (int i = 0; i < n; i++) r(i) = e.vrr[i];
    }
    
    return info;
}

int Eigen::eig(const RealVector &p, const FluxFunction *f, const AccumulationFunction *g, std::vector<eigenpair> &e){
    int n = p.size();

    JetMatrix fm(n), gm(n);
    f->jet(p, fm, 1);
    g->jet(p, gm, 1);

    int info = eig(n, fm.Jacobian().data(), gm.Jacobian().data(), e);

    return info;
}

void Eigen::print_eigen(const vector<eigenpair> &ve){
    for (int i = 0; i < ve.size(); i++){
        printf("Eigenpair #%d\n", i);
        printf("    Eigenvalue = % g + % gi\n", ve[i].r, ve[i].i);
        printf("\n");

        printf("    Left eigenvector =\n");
        for (int j = 0; j < ve[0].vlr.size(); j++){
            printf("    [% g + % gi]\n", ve[i].vlr[j], ve[i].vli[j]);
        }
        printf("\n");

        printf("    Right eigenvector =\n");
        for (int j = 0; j < ve[0].vlr.size(); j++){
            printf("    [% g + % gi]\n", ve[i].vrr[j], ve[i].vri[j]);
        }
        printf("\n");
    }

    return;
}

void Eigen::eps(double e){
    epsilon = e;
    return;
}

double Eigen::eps(void){
    return epsilon;
}

void Eigen::fill_eigenpairs(const FluxFunction *f, const AccumulationFunction *a, const RealVector &u, std::vector<eigenpair> &e){
    // Find the eigenpairs...
    //
    int n = u.size();

    DoubleMatrix JF(n, n);
    DoubleMatrix JG(n, n);
    f->fill_with_jet(n, u.components(), 1, 0, JF.data(), 0);
    a->fill_with_jet(n, u.components(), 1, 0, JG.data(), 0);

    // ...and fill the output.
    //
    eig(n, JF.data(), JG.data(), e);

    return;
}

void Eigen::fill_eigenvalues(const FluxFunction *f, const AccumulationFunction *a, const RealVector &u, std::vector<double> &lambda){
    // Find the eigenpairs...
    //
    std::vector<eigenpair> e;
    fill_eigenpairs(f, a, u, e);

    // ...and fill the output.
    //
    int n = e.size();
    lambda.resize(n);
    for (int i = 0; i < n; i++) lambda[i] = e[i].r;

    return;
}


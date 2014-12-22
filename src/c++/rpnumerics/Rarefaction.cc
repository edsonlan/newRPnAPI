#include "Rarefaction.h"
#include "RarefactionCurve.h"

FluxFunction *Rarefaction::fluxfunction;
AccumulationFunction *Rarefaction::accumulationfunction;
int Rarefaction::type;
int Rarefaction::family;
const Boundary *Rarefaction::boundary;

#ifdef TEST_STENCIL
std::ofstream Rarefaction::stencil;
#endif

double Rarefaction::ddot(int n, double *x, double *y) {
    double p = 0.0;

    for (int i = 0; i < n; i++) p += x[i] * y[i];

    return p;
}

// C = A*B
// A = m times p
// B = p times n
// C = m times n

void Rarefaction::matrixmult(int m, int p, int n, double *A, double *B, double *C) {
    double sum;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            sum = 0.0;
            for (int k = 0; k < p; k++) sum += A[i * p + k] * B[k * n + j];
            C[i * n + j] = sum;
        }
    }

    return;
}

void Rarefaction::fill_with_jet(const RpFunction *flux_object, int n, double *in, int degree, double *F, double *J, double *H) {
    RealVector r(n);
    double *rp = r;
    for (int i = 0; i < n; i++) rp[i] = in[i];

    // Will this work? There is a const somewhere in fluxParams.
    //FluxParams fp(r);
    //flux_object->fluxParams(FluxParams(r)); // flux_object->fluxParams(fp);

    WaveState state_c(r);
    JetMatrix c_jet(n);

    flux_object->jet(state_c, c_jet, degree);

    // Fill F
    if (F != 0) for (int i = 0; i < n; i++) F[i] = c_jet.get(i);

    // Fill J
    if (J != 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                J[i * n + j] = c_jet.get(i, j);
            }
        }
    }

    // Fill H
    if (H != 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    //H[(i * n + j) * n + k] = c_jet(i, j, k); // Check this!!!!!!!!
                    H[i * n + j + n * n * k] = c_jet.get(k, i, j); // This works for the convention adopted in the FluxFunction::jet().
                }
            }
        }
    }

    return;
}

// The flux that is being integrated.
//
// 2012/02/07. ***ELIPTIC REGION***
// The following modification is being considered:
//
// So far, when the rarefaction curve enters the family's eliptic region (which can be determined because the eigenvalue
// associated to the family becomes complex and no longer is a pure real), the flux warns LSODE that an error occurred.
// LSODE then transmits upward the error message, before computing anything else, and halts.
// 
// Now, in order to obtain approximately a point of the rarefaction curve as close as possible to the eliptic region,
// when the flux detects the aforementioned passage from real to complex, it will NOT warn LSODE that there is a problem.
// Instead, it will return to LSODE the real part of the right-eigenvector correspondig to the given family AND will
// modify the first element of param (the one that indicates the family) to reflect this situation, say it can set param[0] = -1.0.
// Then, Rarefaction::curve() must be modified accordingly to check, after each call to LSODE, if
// param[0] == -1.0. If so, along the segment subtended by the previous point (outside the eliptic region) and the new point
// (inside the eliptic region), the zero of the problem's discriminant is to be found by bisection. The point thus found
// will be a close approximation of the intersection between the rarefaction curve and the eliptic region.
//
// NOTE: The discriminant is: det(lambda*DG - DF), where DG is the Jacobian of the accumulation, DF is the Jacobian of the flux.
//

int Rarefaction::flux(int *neq, double *xi, double *in, double *out, int *nparam, double *param) {
    RealVector p_in(*neq, in);

    #ifdef TEST_STENCIL
    stencil << "Flux was called by the solver for point: " << p_in << std::endl;
    #endif

    if (!Rarefaction::boundary->inside(p_in)){
        std::cout << "Rarefaction::flux: The point given " << p_in << " lies OUTSIDE the domain." << std::endl;

        #ifdef TEST_STENCIL
        stencil << "    Error detected: point lies outside the domain!" << std::endl;
        #endif

        return ABORTED_PROCEDURE;
    }

    // The dimension of the problem:
    int n = *neq;

    // The family:
    //    int family = (int)param[0];

    // The reference eigenvector:
    double rev[n];
    for (int i = 0; i < n; i++) rev[i] = param[1 + i];

    int info;
    std::vector<eigenpair> e;

    // Fill the Jacobian
    double FJ[n][n];
    double FH[n][n][n];

    fill_with_jet((RpFunction*) Rarefaction::fluxfunction, n, in, 2, 0, &FJ[0][0], &FH[0][0][0]);

    if (type == RAREFACTION_SIMPLE_ACCUMULATION)
        info = Eigen::eig(n, &FJ[0][0], e);

    else {
        double GJ[n][n];
        double GH[n][n][n];

        fill_with_jet((RpFunction*) Rarefaction::accumulationfunction, n, in, 2, 0, &GJ[0][0], &GH[0][0][0]);
        info = Eigen::eig(n, &FJ[0][0], &GJ[0][0], e);
    }

    if (Rarefaction::family > e.size() - 1){
        std::cout << "Flux. Insufficient eigenpairs! Family = " << Rarefaction::family << "e.size() = " << e.size() << std::endl;
        return ABORTED_PROCEDURE;
    }

    for (int i = 0; i < n; i++) out[i] = e[Rarefaction::family].vrr[i];

    // Check for stop criteria
    // TODO This section is to be tuned.
    if (info == 0) {
        // All eigenvalues must be real. This is true if the number of equations is 2, since
        // one complex eigenvalue implies that the other eigenvalue is complex too.
        // It is necessary to think what to do for higher dimensions. One possibility
        // is to check ONLY if the current family's eigenvalue is complex. According to
        // Dan, this must be thought carefuly.
        for (int i = 0; i < e.size(); i++) {
            if (fabs(e[i].i) > 0) {
                return COMPLEX_EIGENVALUE;

                // ***ELIPTIC REGION***
                // When the eliptic region is reached the eigenvalue is complex.
                // Return the real part of the right eigenvector, as usual, AND
                // make param[0] = -1.0 (or some error code to be defined).
                // DO NOT RETURN COMPLEX_EIGENVALUE, INSTEAD RETURN SUCCESSFUL_PROCEDURE.
                // Thus LSODE will work as usual.
            }
        }

        //        // This family's eigenvalue must be real
        //        if (fabs(e[Rarefaction::family].i) > 0) {
        //            return COMPLEX_EIGENVALUE;
        //        }

        // The following case is unusual, according to Dan, and
        // the action taken must be better sorted out. (2012/02/07)
        //
        // All eigenvalues must be different.
        // This test can be performed thus because the eigencouples are ordered
        // according to the real part of the eigenvalue. Thus, if two eigenvalues
        // are equal, one will come after the other.
        for (int i = 0; i < e.size() - 1; i++) {
            if (e[i].r == e[i + 1].r) {
                return ABORTED_PROCEDURE;
            }
        }

    }
    else {
        return ABORTED_PROCEDURE;
    }

    // The eigenvector to be returned is the one whose inner product with the
    // reference vector is positive.
    if (ddot(n, &(e[Rarefaction::family].vrr[0]), &rev[0]) > 0) {
        for (int i = 0; i < n; i++) out[i] = e[Rarefaction::family].vrr[i];
    }
    else {
        for (int i = 0; i < n; i++) out[i] = -e[Rarefaction::family].vrr[i];
    }

    //        // STOP CRITERION:
    //        // The identity in Proposition 10.11 of
    //        // "An Introduction to Conservation Laws:
    //        // Theory and Applications to Multi-Phase Flow" must not change
    //        // sign (that is, the rarefaction is monotonous).
    //    
    //        double res[n];
    //    
    //        applyH(n, &(e[family].vrr[0]), &H[0][0][0], &(e[family].vrr[0]), &res[0]);
    //        double dlambda_dtk = ddot(n, &res[0], &(e[family].vlr[0])) /
    //                ddot(n, &(e[family].vlr[0]), &(e[family].vrr[0]));
    //        if (dlambda_dtk * ref_speed < 0) {
    //            return ABORTED_PROCEDURE;
    //        }

    #ifdef TEST_STENCIL
    RealVector output(*neq, out);
    stencil << "    No errors detected. Output: " << output << std::endl;
    #endif

    return SUCCESSFUL_PROCEDURE;
}
// Compute the last point of the rarefaction curve when the monotonicity of the eigenvalues
// is violated. This function is totally wrong now.
//
//     previous_point, and new_point
//
// must be of size n + 1. The last component is lambda.
//  

int Rarefaction::compute_last_point(const RealVector &previous_point, const RealVector &new_point, RealVector &last_point) {
    int n = previous_point.size() - 1;

    //last_point.resize(n + 1);

    int info = rar_last_point(n, previous_point, new_point, last_point);
    printf("Inside compute_last_point:\n");

    printf("previous_point = (");
    for (int i = 0; i < n; i++) {
        printf("%g", previous_point.component(i));
        if (i <= (n - 2)) printf(", ");
    }
    printf(")\n");

    printf("new_point = (");
    for (int i = 0; i < n; i++) {
        printf("%g", new_point.component(i));
        if (i <= (n - 2)) printf(", ");
    }
    printf(")\n");

    printf("last_point = (");
    if (info == SUCCESSFUL_PROCEDURE) {
        for (int i = 0; i < n; i++) {
            printf("%g", last_point.component(i));
            if (i <= (n - 2)) printf(", ");
        }
        printf(")\n");
    } else printf("GARBAGE HERE)\n");

    return info;
}

// Compute all the eigenpairs.
//

void Rarefaction::compute_all_eigenpairs(int n, const RealVector &in, std::vector<eigenpair> &e) {
    double p[n];
    for (int i = 0; i < n; i++) p[i] = in.component(i);

    if (type == RAREFACTION_SIMPLE_ACCUMULATION) {
        double FJ[n][n];
        fill_with_jet((RpFunction*) Rarefaction::fluxfunction, n, p, 1, 0, &FJ[0][0], 0);
        Eigen::eig(n, &FJ[0][0], e);
    } else {
        double FJ[n][n], FG[n][n];
        fill_with_jet((RpFunction*) Rarefaction::fluxfunction, n, p, 1, 0, &FJ[0][0], 0);
        fill_with_jet((RpFunction*) Rarefaction::accumulationfunction, n, p, 1, 0, &FG[0][0], 0);
        Eigen::eig(n, &FJ[0][0], &FG[0][0], e);
    }

    return;
}

// Compute the value of the eigenvalue and right eigenvector at a given point for a given family.
// Normally I would use n = in.size(), but this way I can feed this function with a
// vector that contains > n components and even so use its first n components.
//

void Rarefaction::compute_eigenpair(int n, const RealVector &in, double &lambda, RealVector &eigenvector) {

    std::vector<eigenpair> e;
    compute_all_eigenpairs(n, in, e);

    lambda = e[Rarefaction::family].r;

    eigenvector.resize(e[Rarefaction::family].vrr.size());
    for (int i = 0; i < e[Rarefaction::family].vrr.size(); i++) eigenvector.component(i) = e[Rarefaction::family].vrr[i];


    return;
}

// Compute the value of lambda at a given point, a mere wrapper for compute_eigenpair().
// Normally I would use n = in.size(), but this way I can feed this function with a
// vector that contains > n components and even so use its first n components.
//

double Rarefaction::compute_lambda(int n, const RealVector &in) {
    double lambda;

    RealVector eigenvector; // Discarded on exit.

    compute_eigenpair(n, in, lambda, eigenvector);

    return lambda;
}

// Initialize the rarefaction, or, find the second point in the curve.
// The eigenvalue at said point will be stored in the last component of second_point.
//

int Rarefaction::init(const RealVector &initial_point, int increase, double deltaxi, RealVector &second_point) {
    int n = initial_point.size();

    // Eigenvalue and right eigenvector at the initial point.
    double lambda;
    RealVector ev;

    compute_eigenpair(n, initial_point, lambda, ev);

    // Eigenvalues at the candidate points (initial_point +/- deltaxi*ev).
    // Prefixes: m = -, p = +.
    RealVector m(n), p(n);
    for (int i = 0; i < n; i++) {
        m.component(i) = initial_point.component(i) - deltaxi * ev.component(i);
        p.component(i) = initial_point.component(i) + deltaxi * ev.component(i);
    }

    double mlambda = compute_lambda(n, m);
    double plambda = compute_lambda(n, p);

    // Fill the second point of the curve and the eigenvalue thereat.
    second_point.resize(n + 1);

    if (increase == RAREFACTION_SPEED_SHOULD_INCREASE) {
        if (mlambda < lambda && lambda < plambda) {
            for (int i = 0; i < n; i++) second_point.component(i) = p.component(i);
            second_point.component(n) = plambda;
        } else if (mlambda > lambda && lambda > plambda) {
            for (int i = 0; i < n; i++) second_point.component(i) = m.component(i);
            second_point.component(n) = mlambda;
        } else return RAREFACTION_INIT_FAILURE;
    } else if (increase == RAREFACTION_SPEED_SHOULD_DECREASE) {
        if (mlambda < lambda && lambda < plambda) {
            for (int i = 0; i < n; i++) second_point.component(i) = m.component(i);
            second_point.component(n) = mlambda;
        } else if (mlambda > lambda && lambda > plambda) {
            for (int i = 0; i < n; i++) second_point.component(i) = p.component(i);
            second_point.component(n) = plambda;
        } else return RAREFACTION_INIT_FAILURE;
    }

    return RAREFACTION_INIT_OK;
}

// Adapted from the FORTRAN original (at eigen.F, code provided by Dan).
//
// Adapted to the GENERALIZED case by Helmut.
//
// This function computes the directional derivatives for
// generalized problem.
//
//          n: Dimension of the space.
//          p: Point where the directional derivative is computed.
//
// The function returns the directional derivative at p, for
// a given flux and accumulation.
//
// Let l and r be the left- and right-eigenvectors of the GENERALIZED
// SYSTEM OF CONSERVATION LAWS. (notice the difference between this and
// the Stone case), at point p. Let lambda be the associated eigenvalue
// (All of them corresponding to the same family).
//
// Let B the Jacobian of the (NON-TRIVIAL) accumulation.
// Let A the Jacobian of the (NON-TRIVIAL) flux.
//
// Let H be the Hessian of said flux, and M the Hessian of the accumulation.
// Then the (GENERALIZED) directional derivative is (see Chapter
// Numerical Methods thesis Helmut):
//
//     dirdrv = ddot( l, H(r,r) - lambda M(r,r) )/ddot(l, Br).
//
// In particular, D = r^T*H is computed thus:
//
//     D[k] = r^T*H[k],     k = 0,..., (n - 1),
//
// where D[k] is the k-th row of n*n matrix D and H[k] is the
// k-th "matrix" of the Hessian.  (Where there are n Hessians.)
// i.e. it is the Hessian of the k-th component function.
// r^T denotes the transpose of r.
//
// In particular, E = r^T*M is computed thus:
//
//     E[k] = r^T*M[k],     k = 0,..., (n - 1),
//
// where E[k] is the k-th row of n*n matrix E and M[k] is the
// i-th "matrix" of the Hessian.  (Where there are n Hessians.)
//

//double Rarefaction::dirdrv(int n, const RealVector &p){
//    double point[n];
//    for (int i = 0; i < n; i++) point[i] = p.component(i);
//    double A[n][n];

//    double B[n][n];
//    double H[n][n][n];
//    double M[n][n][n];
//    fill_with_jet((RpFunction*)accumulationfunction, n, point, 2, 0, &B[0][0], &M[0][0][0]);
//    fill_with_jet((RpFunction*)fluxfunction,         n, point, 2, 0, &A[0][0], &H[0][0][0]);

//    // Extract the left and right eigenvalues of the generalized system.
//    std::vector<eigenpair> e;
//    int info = Eigen::eig(n, &A[0][0], &B[0][0], e);

//    // Extract the indx-th left and right-eigenvector of the GENERALIZED
//    // PROBLEM (A - lambda B)r=0  and  l(A - lambda B)=0

//    double l[n], r[n];

//    for (int i = 0; i < n; i++){
//        l[i] = e[family].vlr[i];
//        r[i] = e[family].vrr[i];
//    }

//    // Extract lambda.
//    // The i-th eigenvalue must be real. 
//    // The eigenvalues must be chosen carefully in the n-dimensional case.
//    // ALL eigenvalues must be real. Extend this by using a for cycle.
//    //
//    double lambda;

//    if (e[family].i != 0){
//        printf("Inside dirdrv(): Init step, eigenvalue %d is complex: % f %+f.\n", family, e[family].r, e[family].i);
//        return ABORTED_PROCEDURE;     
//    }
//    else lambda = e[family].r;
//    // Compute D and E
//    double D[n][n];
//    double E[n][n];
//    for (int k = 0; k < n; k++){
//       double SubH[n][n];
//       double SubM[n][n];

//       for (int i = 0; i < n; i++){
//           for (int j = 0; j < n; j++){

//               SubH[i][j] = H[k][i][j];
//               SubM[i][j] = M[k][i][j];
//           }
//       }

//       double rowh[n];
//       double rowm[n];

//       matrixmult(1, n, n, &r[0], &SubH[0][0], &rowh[0]);
//       matrixmult(1, n, n, &r[0], &SubM[0][0], &rowm[0]);

///*     matrixmult(n, n, 1, &SubH[0][0], &r[0], &rowh[n]);
//       matrixmult(n, n, 1, &SubK[0][0], &r[0], &rowk[n]); */

//       for (int j = 0; j < n; j++) D[k][j] = rowh[j];             
//       for (int j = 0; j < n; j++) E[k][j] = rowm[j];         
//    }

//    // Compute D*r and E*r
//    double Dtimesr[n];

//    double Etimesr[n];compute_lambda(int n, const RealVector &in)
//    matrixmult(n, n, 1, &D[0][0], &r[0], &Dtimesr[0]);
//    matrixmult(n, n, 1, &E[0][0], &r[0], &Etimesr[0]);

//    // Compute B*r
//    double Btimesr[n];
//    matrixmult(n, n, 1, &B[0][0], &r[0], &Btimesr[0]);

//    // Result
//    double Dtimesr_minus_lambdaEtimesr[n];
//    for (int i = 0; i < n; i++) Dtimesr_minus_lambdaEtimesr[i] = Dtimesr[i] - lambda*Etimesr[i];

//    //for (int i = 0; i < n; i++) printf("l[%d] = %g, r[%d] = %g, Dtimesr[%d] = %g, lambda = %g, Etimesr[%d] = %g\n", i, l[i], i, r[i], i, Dtimesr[i], lambda, i, Etimesr[i]);

//    //printf("Dirdrv = %g divided by %g.\n", ddot(n, Dtimesr_minus_lambdaEtimesr,l), ddot(n, Btimesr, l));

//    return ddot(n, Dtimesr_minus_lambdaEtimesr, l)/ddot(n, Btimesr, l);
////   return ddot(n, Dtimesr, l)/ddot(n, r, l);
//}

double Rarefaction::dirdrv(int n, const RealVector &p, const RealVector &direction) {
    double point[n], dir[n];

//    std::cout << "Dirdrv. Before. n = " << n << std::endl;

    for (int i = 0; i < n; i++) {
        point[i] = p.component(i);
        dir[i] = direction.component(i);
    }

//    std::cout << "Dirdrv. After conversion." << std::endl;

//    cout <<"Ponto: "<<p<<" direcao: "<<direction<<endl;
    return dirdrv(n, &point[0], &dir[0]);
}

// TODO: This method's signature must change. The return value must be an integer
// to let the user know if an error was detected when finding the eigenpairs, or if the point lies outside the domain (etc). Therefore, the dirdrv must be returned as a parameter.
//
// TODO: Remember that this dirdrv is only valid when the problem is d_dt(G) + d_dx(F) = 0.
// For the problem d_dt(G) + d_dx(u*F) = 0 a different dirdrv is needed. Rarefaction will become an abstract class and derived classes will be written to accomodate this situation.

double Rarefaction::dirdrv(int n, double *point, double *dir) {
    int fam = family;

    double A[n][n];
    double B[n][n];

    double H[n][n][n];
    double M[n][n][n];
    fill_with_jet((RpFunction*) accumulationfunction, n, point, 2, 0, &B[0][0], &M[0][0][0]);
    fill_with_jet((RpFunction*) fluxfunction, n, point, 2, 0, &A[0][0], &H[0][0][0]);

    // Extract the left and right eigenvalues of the generalized system.
    std::vector<eigenpair> e;
    int info = Eigen::eig(n, &A[0][0], &B[0][0], e);
//    std::cout << "Dirdrv. e.size() = " << e.size() << std::endl;

    // Extract the indx-th left and right-eigenvector of the GENERALIZED
    // PROBLEM (A - lambda B)r=0  and  l(A - lambda B)=0

    double l[n], r[n];
    double norm = 0.0;
    double dirdrv = 0.0;

    for (int i = 0; i < n; i++) {
        l[i] = e[fam].vlr[i];
        r[i] = e[fam].vrr[i];
    }

    if (ddot(n, r, dir) < 0.0) {
        for (int i = 0; i < n; i++) r[i] = -r[i];
    }

    // Extract lambda.
    // The i-th eigenvalue must be real. 
    // The eigenvalues must be chosen carefully in the n-dimensional case.
    // ALL eigenvalues must be real. Extend this by using a for cycle.
    //
    double lambda;

    if (e[fam].i != 0) {
        printf("Inside dirdrv(): Init step, eigenvalue %d is complex: % f %+f.\n", fam, e[fam].r, e[fam].i);
        return ABORTED_PROCEDURE;
    } else lambda = e[fam].r;

    // TODO: Extract these matrices from the JetMatrix's Hessian.
    double SubH[n][n];
    double SubM[n][n];

    // Nested loops, constructing SubH, SubM as H, M times r and the norm.
    for (int k = 0; k < n; k++) {
        for (int m = 0; m < n; m++) {
            SubH[k][m] = 0.0;
            SubM[k][m] = 0.0;
            norm += l[k] * B[k][m] * r[m];
            for (int i = 0; i < n; i++) {
                SubH[k][m] += H[k][m][i] * r[i];
                SubM[k][m] += M[k][m][i] * r[i];
            }
            // For trivial accumulation, the directional derivative may be simplified as
            //                 dirdrv += l[k]*SubH[k][m]*r[m];
            dirdrv += l[k]*(SubH[k][m] * r[m] - lambda * SubM[k][m] * r[m]);
        }
    }

//    std::cout << "Dirdrv. r = " << RealVector(n, r) << ", lambda = " << lambda << std::endl;
//    std::cout << "    l = " << RealVector(n, l) << std::endl;
//    std::cout << "    Non-null in the Hessian: H[0][0][0] = " << H[0][0][0] << std::endl;
//    std::cout << "    dirdrv = " << dirdrv << ", norm = " << norm << std::endl;
//    std::cout << "    What will be outputed: " << dirdrv / norm << std::endl;

    return dirdrv / norm;
}

// Compute the directional derivative at the rarefaction's initial point.
//
// This function replaces function
//
//     init(const RealVector &initial_point, int increase, double deltaxi, RealVector &second_point),
//
// used heretofore, since it is more robust.
// 
// The user decides if the eigenvalues must increase or decrease as the rarefaction is computed.
// From now on it is assumed that eigenvalues are expected to increase.d_dt(G) + d_dx(F) = 0.
// If so, the directional derivative at the initial point is positive.
// The right eigenvector at the initial point for the current family is returned by LAPACK
// with the correct direction, but since x and y = -x are both eigenvectors for the same
// eigenvalue, it is not to be expected that the eigenvector is also correctly oriented.
// (This also applies to the left eigenvector, but is immaterial in this case, since a said
// eigenvector is multiplied twice, thus effectively cancelling the effects of its orientation.)
// Thus, when computing the directional derivative, it can be affirmed that its absolute value
// is correct, but its sign can be wrong. However, since it is assumed that eigenvalues are
// to increase, the directional derivative is to be positive.
// This being the case, it is so because the right eigenvector was returned by LAPACK with
// the appropriately oriented. If, however, the directional derivative is negative, 
// both the directional derivative and the right eigenvector must be negated.
//

int Rarefaction::initial_dirdrv(int n, const RealVector &p, int increase, double &dd, RealVector &dir) {
    int fam = family;
    double point[n];
    for (int i = 0; i < n; i++) point[i] = p.component(i);
    //
    double A[n][n];
    double B[n][n];
    //
    //    double H[n][n][n];
    //    double M[n][n][n];
    fill_with_jet((RpFunction*) accumulationfunction, n, point, 1, 0, &B[0][0], 0);
    fill_with_jet((RpFunction*) fluxfunction, n, point, 1, 0, &A[0][0], 0);
    //
    //    // Extract the left and right eigenvalues of the generalized system.
    std::vector<eigenpair> e;
    int info = Eigen::eig(n, &A[0][0], &B[0][0], e);

    //    // Extract the indx-th left and right-eigenvector of the GENERALIZED
    //    // PROBLEM (A - lambda B)r=0  and  l(A - lambda B)=0
    //
    double l[n], r[n];
    //    double norm   = 0.0;
    //    double dirdrv = 0.0;
    //
    for (int i = 0; i < n; i++) {
        l[i] = e[fam].vlr[i];
        r[i] = e[fam].vrr[i];
    }
    //
    //    // Extract lambda.
    //    // The i-th eigenvalue must be real.
    //    // The eigenvalues must be chosen carefully in the n-dimensional case.
    //    // ALL eigenvalues must be real. Extend this by using a for cycle.
    //    //
    //    double lambda;
    //
    //    if (e[fam].i != 0){
    //        printf("Inside dirdrv(): Init step, eigenvalue %d is complex: % f %+f.\n", fam, e[fam].r, e[fam].i);
    //        return ABORTED_PROCEDURE;
    //    }
    //    else lambda = e[fam].r;
    //
    //    double SubH[n][n];
    //    double SubM[n][n];
    //
    //    // Nested loops, constructing SubH, SubM as H, M times r and the norm.
    //    for (int k = 0; k < n; k++) {
    //        for (int m = 0; m < n; m++) {
    //            SubH[k][m] = 0.0;
    //            SubM[k][m] = 0.0;
    //            norm += l[k] * B[k][m] * r[m];
    //            for (int i = 0; i < n; i++) {
    //                SubH[k][m] += H[k][m][i]*r[i];
    //                SubM[k][m] += M[k][m][n]*r[n];
    //            }
    //            // For trivial accumulation, the directional derivative may be simplified as
    ////                 dirdrv += l[k]*SubH[k][m]*r[m];
    //            dirdrv += l[k]*(SubH[k][m]*r[m] - lambda*SubM[k][m]*r[m]);
    //        }
    //    }
    //
    //    dd = dirdrv/norm;
    //    dir.resize(n);
    //dd = dirdrv(n, point, dir);

    dd = dirdrv(n,point,r);
    
    if (increase == RAREFACTION_SPEED_SHOULD_INCREASE) {
        if (dd > 0.0) for (int i = 0; i < n; i++) dir.component(i) = r[i];
        else {
            for (int i = 0; i < n; i++) dir.component(i) = -r[i];
            dd = -dd;
        }
    } else if (increase == RAREFACTION_SPEED_SHOULD_DECREASE) {
        if (dd < 0.0) for (int i = 0; i < n; i++) dir.component(i) = r[i];
        else {
            for (int i = 0; i < n; i++) dir.component(i) = -r[i];
            dd = -dd;
        }
    }

//    printf("initial_dirdrv. dd = %lf\n", dd);
    std::cout << "Initial dirdrv: " << dd << std::endl;
    std::cout << "Increase: " << increase << std::endl;
    std::cout << "  r = " << RealVector(n, r) << std::endl;    
    std::cout << "dir = " << dir << std::endl;

    return RAREFACTION_INIT_OK;
}

// Function rar_last_point:
//
// Given two points (and the eigenvalues) p0 and p1, this function computes d0 and d1,
// the directional derivatives at said points for a given family.
// 
// If sign(d0) == -sign(d1) then the function proceeds to compute the point
// where the directional derivative is zero: that is, the point
// where an extremum is reached.
// 
// The resulting point, along with its eigenvalue, is stored in out.
//
// The function returns:
//     -1: d0 and d1 share the same sign, thus the extremum was not found (and out contains garbage).
//      0: Successfuly found the point, out is usable.
//
//int Rarefaction::rar_last_point(int n, const RealVector &p0, const RealVector &p1, RealVector &out){
//    double d0 = dirdrv(n, p0);
//    double d1 = dirdrv(n, p1);

//    if (d1*d0 >= 0.0) return ABORTED_PROCEDURE;

//    double alpha = d1/(d1 - d0); printf("Inside rar_last_point(): alpha = %g\n", alpha);

//    out.resize(n + 1);

//    // This way it is possible to compute the eigenvalue for the last point of the curve.
//    // If the eigenvalue is computed using compute_lambda, a segfault occurs, because
//    // the Jacobian is singular.
//    //
//    for (int i = 0; i <= n; i++) out.component(i) = alpha*p0.component(i) + (1.0 - alpha)*p1.component(i);

//    return SUCCESSFUL_PROCEDURE;
//}

int Rarefaction::rar_last_point(int n, const RealVector &p0, const RealVector &p1, RealVector &out) {

    // Assume p1 = new_point, p0 = previous_point.
    //
    RealVector r_direction(n);
    for (int i = 0; i < n; i++) r_direction.component(i) = p1.component(i) - p0.component(i);

    double d0 = dirdrv(n, p0, r_direction);
    double d1 = dirdrv(n, p1, r_direction);

    if (d1 * d0 >= 0.0) return ABORTED_PROCEDURE;

    // Assume d0 > 0 and d1 < 0. If not, swap.
    RealVector ptemp0(n), ptemp1(n);
    if (d0 > 0.0) {
        for (int i = 0; i < n; i++) {
            ptemp0.component(i) = p0.component(i);
            ptemp1.component(i) = p1.component(i);
        }
    } else {
        for (int i = 0; i < n; i++) {
            ptemp0.component(i) = p1.component(i);
            ptemp1.component(i) = p0.component(i);

            double temp = d0;
            d0 = d1;
            d1 = temp;
        }
    }

    RealVector pmean(n);
    double mean_dirdrv;

    double epsilon = 1e-10;
    int it = 0;
    while (fabs(d0 - d1) > epsilon && it < 100) {
        it++;

        for (int i = 0; i < n; i++) pmean.component(i) = .5 * (ptemp0.component(i) + ptemp1.component(i));
        mean_dirdrv = dirdrv(n, pmean, r_direction);

        if (mean_dirdrv >= 0.0) {
            d0 = mean_dirdrv;
            for (int i = 0; i < n; i++) ptemp0.component(i) = pmean.component(i);
        } else {
            d1 = mean_dirdrv;
            for (int i = 0; i < n; i++) ptemp1.component(i) = pmean.component(i);
        }
    }
    // aqui deberiamos cerciorarnos de que el punto de salida este despues de la infleccion y no antes, o sea que el ultimo mean_dirdrv tenga el mismo signo que en el primer punto de la rarefaccion.
    out.resize(n + 1);
    for (int i = 0; i < n; i++) out.component(i) = pmean(i);

    return SUCCESSFUL_PROCEDURE;
}

int Rarefaction::curve(const RealVector &initial_point,
        int initialize,
        const RealVector *initial_direction,
        int curve_family,
        int increase,
        int type_of_rarefaction,
        ODE_Solver *odesolver,
        double deltaxi,
        const FluxFunction *ff, const AccumulationFunction *aa,
        int type_of_accumulation,
        const Boundary *b,
        std::vector<RealVector> &rarcurve,
        std::vector<RealVector> &inflection_points) {



    // Set the static parameters that will be used throughout.
    // TODO: Decide if increase and deltaxi should be static.
    Rarefaction::fluxfunction = (FluxFunction*) ff;
    Rarefaction::accumulationfunction = (AccumulationFunction*) aa;
    Rarefaction::type = type_of_accumulation;
    Rarefaction::family = curve_family;
    Rarefaction::boundary = b;

    cout << "Rarefaction. Initial point: " << initial_point << endl;


//
//    RealVector test0(3);
//
//    test0.component(0) = 0.209846;
//    test0.component(1) = 0.26888;
//    test0.component(2) = 1.16777;
//
//    RealVector test1(3);
//
//    test1.component(0) = 0.209007;
//    test1.component(1) = 0.269424;
//    test1.component(2) = 1.16346;direction
//
//
//    RealVector directionTest(2);
//    directionTest.component(0) = -0.000839152;
//    directionTest.component(1) = 0.00054363;
//    //
//    //    
//    std::cout << "Test. " << test0 << " dd: " << dirdrv(2, test0, directionTest) << endl;
//    std::cout << "Test2. " << test1 << " dd: " << dirdrv(2, test1, directionTest) << endl;
//
//
//    return -7;

    // Space dimension.
    int n = initial_point.size();

    // These vectors will be used everywhere by the Boundary, and to check the monotonicity of the speed.
    RealVector previous_point(n + 1), new_point(n + 1);

    // Same as previous_point.component(n) and new_point.component(n),
    // but it adds to the legibility of the code.
    double previous_lambda, new_lambda;

    // Clean the output...
    rarcurve.clear();
    inflection_points.clear();

    // ...and store the initial point

    for (int i = 0; i < n; i++) new_point.component(i) = initial_point.component(i);
    new_point.component(n) = compute_lambda(n, initial_point);
    rarcurve.push_back(new_point);

    //    // Initialize the rarefaction and store the second point (lambda is added by init()).
    //    if (initialize == RAREFACTION_INITIALIZE){
    //        int init_info = init(initial_point, increase, deltaxi, new_point);Rarefaction::boundary
    //        if (init_info != RAREFACTION_INIT_OK){
    //            printf("Rarefaction::curve(): Initialization failure.\n");
    //            return init_info;
    //        }
    //    }
    //    else {
    //        RealVector tempev(n);
    //        double templambda;direction
    //        compute_eigenpair(n, new_point, templambda, tempev);

    //        double d = 0;
    //        for (int i = 0; i < n; i++) d += tempev.component(i)*initial_direction->component(i);

    //        printf("d = %f\n", d);
    //        printf("Eigenvector = (");
    //        for (int i = 0; i < n; i++){
    //            printf("%g", tempev.component(i));previous_point
    //            if (i < n - 1) printf(", ");
    //        }
    //        printf(")\n");

    //        if (d >= 0.0) for (int i = 0; i < n; i++) new_point.component(i) += deltaxi*tempev.component(i);
    //        else          for (int i = 0; i < n; i++) new_point.component(i) -= deltaxi*tempev.component(i);

    //        new_point.component(n) = compute_lambda(n, new_point);

    //        printf("New point   = (");
    //        for (int i = 0; i < n; i++){
    //            printf("%g", new_point.component(i));
    //            if (i < n - 1) printf(", ");
    //        }
    //        printf(")\n");
    //    }

    //    rarcurve.push_back(new_point);
    //    new_lambda = new_point.component(n);

    // BEGIN Prepare the parameters to be passed to LSODE //
    int ml; // Not used.
    int mu; // Not used.

    // ???
    int nrpd = 4;

    // Is the tolerance the same for all the elements of U (1) or not (2)?
    int itol = 2; // 1: atol scalar; 2: atol array.
    double rtol = 1e-4;
    double atol[n];
    for (int i = 0; i < n; i++) atol[i] = 1e-6;

    // The Jacobian is provided by the user.
    // int mf = 21; 
    // The Jacobian is NOT provided by the user.
    int mf = 22;
    // Lsode uses rwork to perform its computations.
    // lrw is the declared length of rwork
    int lrw;
    if (mf == 10) lrw = 20 + 16 * n;
    else if (mf == 21 || mf == 22) lrw = 22 + 9 * n + n * n;
    else if (mf == 24 || mf == 25) lrw = 22 + 10 * n + (2 * ml + mu) * n;
    double rwork[lrw];

    // Normal computation of values at tout.
    int itask = 1;

    // Set to 1 initially.
    // This is where LSODE's info parameter. Must be set to 1 the first time.
    int istate = 1;
    // No optional inputs
    int iopt = 0;

    // Lsode uses iwork to perform its computations.direction
    // liw is the declared length of iwork
    int liw;
    if (mf == 10) liw = 20;
    else if (mf == 21 || mf == 22 || mf == 24 || mf == 25) liw = 20 + n;
    int iwork[liw];
    // END   Prepare the parameters to be passed to LSODE //

    // The point LSODE uses.
    double p[n];
    for (int i = 0; i < n; i++) p[i] = new_point.component(i);

    // Independent parameter. Not used by flux(), but needed by LSODE.
    double xi = 0.0, new_xi = deltaxi;

    // Reference vector (passed as param).
    int nparam = n + 1;
    double param[nparam];
    param[0] = (int) family;
    //    for (int i = 0; i < n; i++) param[1 + i] = new_point.component(i) - initial_point.component(i);

    //int info = SUCCESSFUL_PROCEDURE;

    // Compute the correct direction of the right-eigenvector at the initial point.
    // This will solve the problem that arises when the initial point is closer to the
    // inflection than delta_xi.

    RealVector r_direction(n);
    double new_dirdrv, previous_dirdrv;

    if (initialize == RAREFACTION_INITIALIZE) initial_dirdrv(n, initial_point, increase, new_dirdrv, r_direction);
    else {
        for (int i = 0; i < n; i++) r_direction.component(i) = initial_direction->component(i);
        int info_initial = initial_dirdrv(n, initial_point, increase, new_dirdrv, r_direction);
        if (info_initial == ABORTED_PROCEDURE) return ABORTED_PROCEDURE;
    }

    for (int i = 0; i < n; i++) {
        param[1 + i] = r_direction.component(i);
        new_point.component(i) = initial_point.component(i);
    }

    // Compute the curve.
    #ifdef TEST_STENCIL
//    system("rm -f stencil.txt");
    stencil.open("stencil.txt", ios::trunc);
    stencil << "Rarefaction initiated for point " << initial_point << ", family = " << family << ", increase = " << increase << ".\n\n" << std::endl;
    #endif

    while (true) {
        // Added 2011-11-30:
        previous_dirdrv = new_dirdrv;

        // Update the previous point & lambda
        for (int i = 0; i <= n; i++) previous_point.component(i) = new_point.component(i);
        //previous_point.component(n) = previous_lambda = new_point.component(n);
        previous_lambda = new_point.component(n);

        // Invoke LSODE.
//        std::cout << "Curve, before lsode_" << std::endl;

        #ifdef TEST_STENCIL
        RealVector temp_point(n);
        for (int i = 0; i < n; i++) temp_point(i) = previous_point(i);

        stencil << "\n\nLSODE will be called now for point " << temp_point << "." << std::endl;
        #endif

        lsode_(&flux, &n, p, &xi, &new_xi, &itol, &rtol, atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, 0, &mf, &nparam, param);

//        // Euler solver BEGIN
//        EulerSolver es(Rarefaction::boundary, 5);
//        RealVector init_point(n, p), final_point(n);

//        istate = EulerSolver::euler_solver(&flux_wrapper, (int*)0, (int*)param, 
//                                           xi,     init_point,  
//                                           new_xi, final_point,
//                                           (int*)(&es), 0);


//        for (int i = 0; i < n; i++) p[i] = final_point(i);
//        // Euler solver END

        #ifdef TEST_STENCIL
        for (int i = 0; i < n; i++) temp_point(i) = p[i];

        stencil << "LSODE returned the point  " << temp_point << "." << std::endl;
        #endif

        if (istate == ABORTED_PROCEDURE){
            #ifdef TEST_STENCIL
            stencil << "LSODE returned info = " << istate << ". Leaving...\n\n" << std::endl;
            stencil.close();
            #endif

            return ABORTED_PROCEDURE;
        }
        //        printf("LSODE: info = %d\n", istate);

        // ***ELIPTIC REGION***
        // 2012/02/07.
        // Rarefaction::flux() will, upon entering the eliptic region, modify the first element of param,
        // and will pass to LSODE the real part of the right eigenvector. LSODE will return a point
        // inside the eliptic region.
        // Along the segment subtended by the previous point and the new point, the zero of the discriminant of
        // the problem's generalized Jacobian will be found, by bisection. Iteration will halt afterwards.
        // DO IT HERE

        // Update new_point.

        for (int i = 0; i < n; i++) new_point.component(i) = p[i];
        new_point.component(n) = new_lambda = compute_lambda(n, new_point);


        // BEGIN Check Boundary //
        // Modified RectBoundary so that the intersection can be tested using RealVectors of size
        // greater than the dimension of the space the RectBoundary is in.
        int where_out;
        RealVector r;
        int intersection_info = boundary->intersection(previous_point, new_point, r, where_out);

        //        printf("Inside while. previous_point = (");
        //        for (int i = 0; i < n; i++){
        //            printf("%g", previous_point.component(i));
        //            if (i < n - 1) printf(", ");
        //        }
        //        printf(")\n");
        //        printf("Inside while.      new_point = (");
        //        for (int i = 0; i < n; i++){
        //            printf("%g", new_point.component(i));
        //            if (i < n - 1) printf(", ");
        //        }
        //        printf(")\n");

        if (intersection_info == 1) {
            // Both points inside. Carry on with the rest of the tests, etc.
        } else if (intersection_info == 0) {
            // One point is inside, the other is outside. 
            // Store the point lying in the domain's border and get out.
            r.resize(n + 1);
            r.component(n) = compute_lambda(n, r);
            rarcurve.push_back(r);

            printf("Reached boundary\n");

            #ifdef TEST_STENCIL
            stencil << "Reached boundary. Leaving...\n\n" << std::endl;
            stencil.close();
            #endif

            return SUCCESSFUL_PROCEDURE;
        } else {
            // Both points lie outside the domain. Something went awfully wrong here.
            printf("Both outside\n");
            printf("previous_point = (");
            for (int i = 0; i < n; i++) {
                printf("%g", previous_point.component(i));
                if (i < n - 1) printf(", ");
            }
            printf(")\n");

            printf("new_point      = (");
            for (int i = 0; i < n; i++) {
                printf("%g", new_point.component(i));
                if (i < n - 1) printf(", ");
            }
            printf(")\n");

            return ABORTED_PROCEDURE;
        }
        // END   Check Boundary //

        // BEGIN Check for monotonicity //
        //        if (increase != RAREFACTION_SPEED_NEUTRAL){ABORTED_PROCEDURE
        //            if ((new_lambda > previous_lambda && increase == RAREFACTION_SPEED_SHOULD_DECREASE) || 
        //                (new_lambda < previous_lambda && increase == RAREFACTION_SPEED_SHOULD_INCREASE)){
        for (int i = 0; i < n; i++) r_direction.component(i) = new_point.component(i) - previous_point.component(i);

        new_dirdrv = dirdrv(n, new_point, r_direction);
//        std::cout << "Curve, dirdrv computed." << std::endl;

        //            printf("new_dirdrv = %lg, previous_dirdrv = %g, new_dirdrv*previous_dirdrv = %g\n", new_dirdrv, previous_dirdrv, new_dirdrv*previous_dirdrv);
        //            printf("new_lambda = %lg, previous_lambda = %g\n", new_lambda, previous_lambda);

//        std::cout << "Point = " << new_point << std::endl;
//        std::cout << "    prev dirdrv = " << previous_dirdrv << ", new dirdrv = " << new_dirdrv << std::endl;
        if (new_dirdrv * previous_dirdrv <= 0.0) {
            std::cout << "Inflection seems to be near. Will compute bisection."  << std::endl;
            // The code commented below has been superseded by a new approach. Eliminate it after all tests are satisfactory.

//            printf("Ok");
//            // printf("new_lambda = %g; previous_lambda = %g.\n", new_lambda, previous_lambda);

//            // Find the point where lambda reaches a minimum, store it and get out.
//            RealVector last_point;
//            int info_compute_last_point = compute_last_point(previous_point, new_point, last_point);
//            if (info_compute_last_point == SUCCESSFUL_PROCEDURE) {
//                // The inflection will only be added to the rarefaction curve when it is not used
//                // for the integral curve. In that case is better to ommit that point
//                // in order to avoid arrow clutter when displaying the results.
//                // The value of lambda at the inflection point is not being calculated
//                // and this situation affects the method by which the arrows are created.


//                // Pablo mudou essa linha para que não aparecesse uma seta a mais de direcao contrária na rarefacao (Conferir com o Rodrigo)
//                //if (type_of_rarefaction == RAREFACTION_FOR_ITSELF) rarcurve.push_back(last_point);


//                if (type_of_rarefaction == RAREFACTION_FOR_ITSELF) rarcurve.push_back(new_point);


//            }
//            //                else printf("Last point discarded.\n");

            double t_in  = xi;
            double t_fin = xi + deltaxi;

            // This will be changed in the near future. Right now new_point has, as its last component, the eigenvalue.
            // Dan and I have discussed this matter and we both agree that the lambdas should be moved to a separate array.
            // When that happens, dirdrv will have no need of n, since all the components of the point will be really components
            // and not some extra info on the point.
            //
            RealVector p_in  = previous_point;
            p_in.resize(n);

            RealVector p_fin = new_point;
            p_fin.resize(n);

            double bisection_epsilon = 1e-10;

            // Output here:
            double c_t;
            RealVector p_c;

//            // ODE Solver:
//            EulerSolver es(5);

            int info_bisection = Bisection::bisection_method(t_in,  p_in,
                                                             t_fin, p_fin,
                                                             bisection_epsilon, 
                                                             c_t, p_c,
                                                             //&flux_wrapper, (int*)0, (int*)param,
                                                             &flux, (int*)0, (double*)param,
                                                             odesolver,
                                                             /*&EulerSolver::euler_solver, (int*)(&es), 0,*/
                                                             &inflection_signal_event, 0 /*int *signal_event_object*/, (int*)&r_direction /*int *signal_event_data*/);

            std::cout << "Bisection done. p_c = " << p_c << std::endl;

            // Add the point obtained, notwithstanding the value of info.
            //
            double lambda_p_c = compute_lambda(n, p_c);
            p_c.resize(n + 1);
            p_c(n) = lambda_p_c;

            rarcurve.push_back(p_c);

            if (info_bisection == BISECTION_FUNCTION_ERROR){
                #ifdef TEST_STENCIL
                stencil << "An error was reported by the signal function when called by Bisection. Leaving..." << std::endl;
                stencil.close();
                #endif

                std::cout << "An error was reported by the signal function when called by Bisection. Leaving..." << std::endl;

                return ABORTED_PROCEDURE;
            }

            if (info_bisection == BISECTION_EQUAL_SIGN){
                #ifdef TEST_STENCIL
                stencil << "Bisection detected that the signal event function has the same sign in both points. Leaving..." << std::endl;
                stencil.close();
                #endif

                std::cout << "Bisection detected that the signal event function has the same sign in both points. Leaving..." << std::endl;

                return ABORTED_PROCEDURE;
            }

            if (info_bisection == BISECTION_CONVERGENCE_ERROR){
                #ifdef TEST_STENCIL
                stencil << "Bisection did not converge when computing the inflection. Leaving..." << std::endl;
                stencil.close();
                #endif

                std::cout << "Bisection did not converge when computing the inflection. Leaving..." << std::endl;

                return ABORTED_PROCEDURE;
            }
            else {
                #ifdef TEST_STENCIL
                stencil << "Bisection converged when computing the inflection. Leaving..." << std::endl;
                stencil.close();
                #endif

                std::cout << "Bisection converged when computing the inflection. Leaving..." << std::endl;

                return RAREFACTION_REACHED_INFLECTION;
            }

            //std::cout << "Rarefaction. Inflection point at: " << last_point << std::endl;

            #ifdef TEST_STENCIL
            stencil << "Rarefaction is no longer monotonous:" << std::endl << std::endl;
            stencil << "    Dirdrv was " << previous_dirdrv << ", now is " << new_dirdrv << std::endl << std::endl;
            stencil.close();
            #endif

            printf("RAREFACTION_NOT_MONOTONOUS\n");

            if (type_of_rarefaction == RAREFACTION_FOR_ITSELF) return RAREFACTION_NOT_MONOTONOUS;
            else if (type_of_rarefaction == RAREFACTION_AS_ENGINE_FOR_INTEGRAL_CURVE) {
                // Update the list of inflection points found so far.
                //
//                inflection_points.push_back(last_point);
//                inflection_points.push_back(p_c);
            }
        }

        // Store the point and the eigenvalue and continue.
        //printf("Rarefaction, size = %d\n", rarcurve.size());
        rarcurve.push_back(new_point);
        //        }
        //        else rarcurve.push_back(new_point);
        // END   Check for monotonicity //

        // Update the independent parameters.
        xi = new_xi;
        new_xi += deltaxi;

        // Update the reference vector.
        for (int i = 0; i < n; i++) param[1 + i] = new_point.component(i) - previous_point.component(i);
//        std::cout << "Curve, updated." << std::endl;
    }

    #ifdef TEST_STENCIL
    stencil.close();
    #endif

    //return SUCCESSFUL_PROCEDURE;
}

int Rarefaction::flux_wrapper(const RealVector &in, RealVector &out, int *flux_object, int *flux_data){
    int neq = in.size();
    double xi; // Not used.

    out.resize(neq);

    int nparam = neq + 1;

    int info = flux(&neq, &xi, (double*)in.components(), out.components(), &nparam, (double*)flux_data);

    if (info == SUCCESSFUL_PROCEDURE) return BISECTION_FUNCTION_OK;
    else                              return BISECTION_FUNCTION_ERROR;
}

//int Rarefaction::flux_wrapper_for_lsode(){
//    
//}

//int Rarefaction::lsode_solver(int (*field)(const RealVector &, RealVector &, int *, int *), int *function_object, int *function_data, 
//                              const double init_time,  const RealVector &init_point,  
//                              const double final_time,       RealVector &final_point,
//                              int *lsode_object, int *lsode_data){
//    //
//    int n = init_point.size();

//    int ml; // Not used.
//    int mu; // Not used.

//    // ???
//    int nrpd = 4;

//    // Is the tolerance the same for all the elements of U (1) or not (2)?
//    int itol = 2; // 1: atol scalar; 2: atol array.
//    double rtol = 1e-4;
//    double atol[n];
//    for (int i = 0; i < n; i++) atol[i] = 1e-6;

//    // The Jacobian is provided by the user.
//    // int mf = 21; 
//    // The Jacobian is NOT provided by the user.
//    int mf = 22;
//    // Lsode uses rwork to perform its computations.
//    // lrw is the declared length of rwork
//    int lrw;
//    if (mf == 10)                  lrw = 20 + 16 * n;
//    else if (mf == 21 || mf == 22) lrw = 22 + 9 * n + n * n;
//    else if (mf == 24 || mf == 25) lrw = 22 + 10 * n + (2 * ml + mu) * n;
//    double *rwork = new double[lrw];

//    // Normal computation of values at tout.
//    int itask = 1;

//    // Set to 1 initially.
//    // This is where LSODE's info parameter. Must be set to 1 the first time.
//    int istate = 1;
//    // No optional inputs
//    int iopt = 0;

//    // Lsode uses iwork to perform its computations.
//    // liw is the declared length of iwork
//    int liw;
//    if (mf == 10)                                          liw = 20;
//    else if (mf == 21 || mf == 22 || mf == 24 || mf == 25) liw = 20 + n;
//    int *iwork = new int[liw];
//    // END   Prepare the parameters to be passed to LSODE //

//    // The point LSODE uses.
//    final_point = init_point;

//    // Independent parameter. Not used by flux(), but needed by LSODE.
//    double xi = init_time, new_xi = final_time;

//    // Reference vector (passed as param).
//    int nparam = n + 1;
//    double param[nparam];
//    param[0] = (int) family;

//    lsode_(field, &n, final_point.components(), &xi, &new_xi, &itol, &rtol, atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, 0, &mf, &nparam, param);
//    
//    // Clean up
//    //
//    delete [] iwork;
//    delete [] rwork;

//    if (istate == ABORTED_PROCEDURE) return BISECTION_SOLVER_ERROR;
//    else                             return BISECTION_SOLVER_OK;
//}

int Rarefaction::inflection_signal_event(const RealVector & where, double & directional_derivative_measure, int *signal_object, int *reference_direction){
    // TODO: When the rarefaction is rewritten as non-static, the line below will be rewritten thus (or whereabouts):
//    RealVector *direction = ((Rarefaction*)signal_event_object)->(RealVector*)signal_event_data;

    RealVector *direction = (RealVector*)reference_direction;

    // TODO: Remember that dirdrv's signature will change in the future! This is why an intermediate step
    //       was used here.
    //
//    std::cout << "Inflection_signal, before dirdrv." << std::endl;

    double d = dirdrv(where.size(), where, *direction);

//    std::cout << "Inflection_signal, after dirdrv." << std::endl;

    directional_derivative_measure = d;

//    std::cout << "where.size() = " << where.size() << std::endl;

    // TODO: When dirdrv's signature changes, this method's return value needs not be fixed.
    return BISECTION_FUNCTION_OK;
}


/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) HyperCube.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "HyperCube.h"
/*
 * ---------------------------------------------------------------
 * Definitions:
 */

/*
 * If n is out of range, factorial returns 0.
 */
int HyperCube::factorial(int n) {
    if (n < 0 || n > 7) {
        printf("Error when computing the factorial. n = %d. Will abort now.\n", n);
        return ERROR;
    }
    if (n == 0) return 1;
    else return n * factorial(n - 1);
}

int HyperCube::pow2(int n) {
    int temp = 1;
    int i;

    if (n < 0) return ERROR;

    for (i = 0; i < n; i++) temp *= 2;

    return temp;
}

void HyperCube::Eval_dimf(int *dimFACE) {
    // This table contains the number of distinct m-dimensional faces for all simplices of n-dimensions.
    int i, j;
    int n = 5;
    int m = 6; // Increasing the matrix means increasing these values.

    /* The values for n and m are:
        n \ m |   0  |    1  |     2  |     3  |     4  |     5
       -------+------+-------+--------+--------+--------+--------
          2   |   4  |    5  |     2  |     -  |     -  |     -
          3   |   8  |   19  |    18  |     6  |     -  |     -
          4   |  16  |   65  |   110  |    84  |    24  |     -
          5   |  32  |  211  |   570  |   750  |   480  |   120
          6   |  64  |  665  |  2702  |  5460  |  5880  |  3240
     */

    //initializing the matrix as zero, for following errors.
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            dimFACE[i * m + j] = ERROR;
        }
    }

    dimFACE[0 * m + 0] = 4;
    dimFACE[0 * m + 1] = 5;
    dimFACE[0 * m + 2] = 2;

    dimFACE[1 * m + 0] = 8;
    dimFACE[1 * m + 1] = 19;
    dimFACE[1 * m + 2] = 18;
    dimFACE[1 * m + 3] = 6;

    dimFACE[2 * m + 0] = 16;
    dimFACE[2 * m + 1] = 65;
    dimFACE[2 * m + 2] = 110;
    dimFACE[2 * m + 3] = 84;
    dimFACE[2 * m + 4] = 24;

    dimFACE[3 * m + 0] = 32;
    dimFACE[3 * m + 1] = 211;
    dimFACE[3 * m + 2] = 570;
    dimFACE[3 * m + 3] = 750;
    dimFACE[3 * m + 4] = 480;
    dimFACE[3 * m + 5] = 120;

    dimFACE[4 * m + 0] = 64;
    dimFACE[4 * m + 1] = 665;
    dimFACE[4 * m + 2] = 2702;
    dimFACE[4 * m + 3] = 5460;
    dimFACE[4 * m + 4] = 5880;
    dimFACE[4 * m + 5] = 3240;

    return;
}

//int HyperCube::CubeFace(int dim, int faceDim /*Provalmente faltam ponteros*/) {
//    //In case of ERROR, CubeFace returns the value 0
//    int i, j;
//
//    //inicializing variables
//    int n_ = dim;
//    int m_ = faceDim;
//    int ncvert_ = pow2(n_);
//    int nsimp_ = factorial(n_);
//
//    if (n_ < m_) {
//        printf("Error simplex dimension exceeds geometric dimension.\n");
//        return ERROR;
//    }
//
//    if (n_ < 2 || m_ < 0) {
//        printf("Error, simplex dimension is not valid.\n");
//        return ERROR;
//    }
//    if (n_ > 6 || m_ > 5) {
//        printf("Error, the matrix needs to be bigger (n = %d, m = %d).\n", n_, m_);
//        return ERROR;
//    }
//
//    if (nsimp_ == 0) {
//        printf("Error, dimension exceeds range\n");
//        return ERROR;
//    }
//    //    int dimf_ = dimFACE(n_, m_); // dimFACE(n_); // TODO Falta configurar dimFACE
//    int numberOfCombinations = combination(n_ + 1, m_ + 1);
//
//    //inicializing arrays
//    int storn_[n_ + 1];
//    int storm_[m_ + 1];
//    double cvert_[ncvert_][n_];
//    int bsvert_[n_ + 1][n_];
//    int perm_[n_][nsimp_];
//    int comb_[numberOfCombinations][m_ + 1];
//
//    //inicializing arrays dimensions
//    int nsface_ = mkcomb(&comb_[0][0], n_ + 1, m_ + 1);
//    printf("nsface vale: %d\n", nsface_);
//
//    int fnbr_[nsface_][nsface_];
//
//    int dimFACE[5][6]; // This is for the maximal values n_ = 6 and m_ = 5, see Eval_dimf.
//    Eval_dimf(&dimFACE[0][0]);
//    //TODO: Ver se nface_ eh igual a dimf_ + 1?
//
//    int dimf_ = dimFACE[n_ - 2][m_]; // The matrix dimFACE ranges for
//
//    if (dimf_ == ERROR) {
//        printf("Error, simplices dimensions out of range.\n");
//        return ERROR;
//    }
//
//    //inicializing another arrays, it were globally defined in java
//    int facptr_[nsimp_][nsface_];
//    int face_[m_ + 1][dimf_];
//
//    int np, mp;
//    for (mp = 0; mp < m_ + 1; mp++) {
//        for (np = 0; np < dimf_; np++) {
//            face_[mp][np] = -1;
//        }
//    }
//
//    // funcoes de inicializacao
//    mkcube(&cvert_[0][0], &bsvert_[0][0], &perm_[0][0], ncvert_, nsimp_, n_); //TODO Talvez deve passar comb[][]
//    //Debug!!!!!!!
//
//    /*** In java mkface was a subroutine, here is between the code, if need it it is commented latter. ***/
//    int nface_ = mkface(&face_[0][0], &facptr_[0][0], &fnbr_[0][0], dimf_, nsimp_, n_, m_, nsface_,
//            &bsvert_[0][0], &comb_[0][0], &perm_[0][0], &storn_[0], &storm_[0]);
//    printf("nface vale: %d e dimf vale: %d\n", nface_, dimf_);
//
//
//    //compute the neighbor array
//    //Inicializacao da matriz fnbr_
//    //global variable initialization
//    /*    int fnbr_[nsface_][nsface_];
//        //zero out  fnbr_
//        for (j = 0; j < nsface_; j++) {
//            for (i = 0; i < nsface_; i++) {
//                fnbr_[j][i] = 0;
//            }
//        }
//
//        mkfnbr(&fnbr_[0][0], &comb_[0][0], n_, m_, nsface_, &storn_[0]);
//
//        int nface_ = mkfcfp(&face_[0][0], &facptr_[0][0],
//                        nsface_, dimf_, nsimp_, m_, n_,
//                        &bsvert_[0][0], &comb_[0][0], &perm_[0][0], &storm_[0]);
//        if (nface_ == ERROR) {
//            printf("Error, dimensions do not match. Increase dimf_.\n");
//            return ERROR;
//        }
//     */
//    /*** Here the mkface ends, but needs ***/
//
//
//    //    printf("nface_ = %d, dimf_ = %d\n", nface_, dimf_);
//
//    /*    int ni, nj;
//        printf("face_[%d][%d] = \n", m_+1, dimf_);
//        for (ni = 0; ni < m_ + 1; ni++) {
//            for (nj = 0; nj < dimf_; nj++) {
//                printf("%d ",face_[ni][nj]);
//            }
//            printf("\n");
//        }
//     */
//
//    //Fim do debug
//    //Global variable initialization TODO: NOTAR QUE ESTAVAM EMBAIXO
//    int ptrf_[nface_][2];
//    int ptfatt_[nface_][2];
//
//    int fdim_, faceptInd_;
//    findInds(&faceptInd_, &fdim_, n_, m_, &cvert_[0][0], &face_[0][0], dimf_, nface_, nsimp_, nsface_, &facptr_[0][0]);
//    int facatt_[fdim_][2];
//    int facept_[faceptInd_][2];
//
//
//
//    mkflst(&facept_[0][0], &ptrf_[0][0], &facatt_[0][0], &ptfatt_[0][0],
//            dimf_, fdim_, n_, m_, nface_, nsface_, nsimp_, ncvert_, faceptInd_,
//            &cvert_[0][0], &face_[0][0], &facptr_[0][0], &storm_[0]);
//    //Debug !!!
//
//    //Fim do debug
//
//    //    putmf("CVERT", &cvert_[0][0], ncvert_, n_);
//    //    putmi("BSVERT", &bsvert_[0][0], n_ + 1, n_);
//    //    putmi("PERM", &perm_[0][0], n_, nsimp_);
//    //    putmi("COMB", &comb_[0][0], numberOfCombinations, m_ + 1);
//    //    putmi("FACE", &face_[0][0], m_ + 1, dimf_);
//    //    putmi("FACPTR", &facptr_[0][0], nsimp_, nsface_);
//
//
//    double foncub[m_][ncvert_];
//    for (i = 0; i < ncvert_; i++) {
//        for (j = 0; j < m_; j++) foncub[j][i] = 0;
//    }
//
//    if (n_ == 4 && m_ == 3) {
//        foncub[0][0] = -5.7000e+00;
//        foncub[1][0] = -1.9000e+00;
//        foncub[2][0] = -2.7800e+00;
//
//        foncub[0][1] = -7.0000e-01;
//        foncub[1][1] = -9.0000e-01;
//        foncub[2][1] = -2.6800e+00;
//
//        foncub[0][2] = -2.7000e+00;
//        foncub[1][2] = -9.0000e-01;
//        foncub[2][2] = 2.2200e+00;
//
//        foncub[0][3] = 2.3000e+00;
//        foncub[1][3] = 1.0000e-01;
//        foncub[2][3] = 2.3200e+00;
//
//        foncub[0][4] = -3.7000e+00;
//        foncub[1][4] = 1.0000e-01;
//        foncub[2][4] = -7.8000e-01;
//
//        foncub[0][5] = 1.3000e+00;
//        foncub[1][5] = 1.1000e+00;
//        foncub[2][5] = -6.8000e-01;
//
//        foncub[0][6] = -7.0000e-01;
//        foncub[1][6] = 1.1000e+00;
//        foncub[2][6] = 4.2200e+00;
//
//        foncub[0][7] = 4.3000e+00;
//        foncub[1][7] = 2.1000e+00;
//        foncub[2][7] = 4.3200e+00;
//
//        foncub[0][8] = -4.7000e+00;
//        foncub[1][8] = 1.1000e+00;
//        foncub[2][8] = 2.2000e-01;
//
//        foncub[0][9] = 3.0000e-01;
//        foncub[1][9] = 2.1000e+00;
//        foncub[2][9] = 3.2000e-01;
//
//        foncub[0][10] = -1.7000e+00;
//        foncub[1][10] = 2.1000e+00;
//        foncub[2][10] = 5.2200e+00;
//
//        foncub[0][11] = 3.3000e+00;
//        foncub[1][11] = 3.1000e+00;
//        foncub[2][11] = 5.3200e+00;
//
//        foncub[0][12] = -2.7000e+00;
//        foncub[1][12] = 3.1000e+00;
//        foncub[2][12] = 2.2200e+00;
//
//        foncub[0][13] = 2.3000e+00;
//        foncub[1][13] = 4.1000e+00;
//        foncub[2][13] = 2.3200e+00;
//
//        foncub[0][14] = 3.0000e-01;
//        foncub[1][14] = 4.1000e+00;
//        foncub[2][14] = 7.2200e+00;
//
//        foncub[0][15] = 5.3000e+00;
//        foncub[1][15] = 5.1000e+00;
//        foncub[2][15] = 7.3200e+00;
//    }
//
//
//    int dims_ = 3; //ANTES ESTAVA EM 100!!
//    int dime_ = 20;
//    int nsoln_ = -1;
//    double sol_[n_][dims_];
//    int sptr[nface_];
//    int solptr[nsimp_][nsface_];
//    for (i = 0; i < nsimp_; i++) {
//        for (j = 0; j < nsface_; j++) solptr[i][j] = 0;
//    } //TODO: Revisar como "solptr" eh modificada, os numero sao muito estranhos
//
//    int edges_[2][dime_];
//    for (i = 0; i < 2; i++) {
//        for (j = 0; j < dime_; j++) edges_[i][j] = -6;
//    } //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos
//
//    int smpedg_[nsimp_][2];
//    for (i = 0; i < nsimp_; i++) {
//        for (j = 0; j < 2; j++) smpedg_[i][j] = 0;
//    } //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos
//
//    double u[n_][m_ + 1];
//    double g[m_][m_ + 1];
//    double gx[m_];
//    double x[m_];
//    int wrki[m_];
//    int exstfc[nface_];
//    for (i = 0; i < nface_; i++) exstfc[i] = 1;
//
//    nsoln_ = cubsol(&solptr[0][0], &sol_[0][0], dims_, &sptr[0], nsoln_, &foncub[0][0], &exstfc[0],
//            &face_[0][0], &facptr_[0][0], dimf_, &cvert_[0][0], ncvert_, n_, m_, nsimp_, nsface_, nface_, &u[0][0],
//            &g[0][0], &x[0], &wrki[0]);
//
//    //    putmf("g", &g[0][0], m_, m_ + 1);
//    //
//    //    putmi("FNBR", &fnbr_[0][0], nsface_, nsface_);
//    //    putmf("FONCUB", &foncub[0][0], m_, ncvert_);
//    //    putmf("SOL", &sol_[0][0], n_, dims_);
//    //    putmi("SPTR", &sptr[0], nface_, 1);
//    //    putmi("SOLPTR", &solptr[0][0], nsimp_, nsface_);
//
//    int nedges_;
//    nedges_ = mkedge(&edges_[0][0], dime_, nedges_, &smpedg_[0][0], &solptr[0][0], &fnbr_[0][0], nsimp_, nsface_);
//
//    //    putmi2("EDGES", &edges_[0][0], 2, dime_, nedges_);
//    //    putmi("SMPEDG", &smpedg_[0][0], nsimp_, 2);
//    //
//
//    return OK;
//}

/**
 * @param altFace The new face array
 * @throws IllegalArgumentException If the array has an invalid size.
 */

/*  TODO Parece ser que nao temos como fazer essa mensagem de erro funcionar.
void setFaceArray(int[] [] altFace) {
    //Check sizes

/**
 * This routine creates the standard hypercube vertices (in cvert_),
 * the basic simplex vertices (in bsvert_), and the permutations of 1, 2, ..., n_  (in permutations).
 */

void HyperCube::cpp_mkcube(Matrix<double> &cpp_cvert_, // cvert_[ncvert_][n_]
        Matrix<int> &cpp_bsvert_, // bsvert_[n_ + 1][n_]
        Matrix<int> &cpp_perm_, // perm_[n_][nsimp_]
        int ncvert_, int nsimp_, int n_) {

    double cvert_[ncvert_][n_];
    int bsvert_[n_ + 1][n_];
    int perm_[n_][nsimp_];

    mkcube(&cvert_[0][0], &bsvert_[0][0], &perm_[0][0], ncvert_, nsimp_, n_);

    cpp_cvert_ = Matrix<double>(ncvert_, n_, &cvert_[0][0]);
    cpp_bsvert_ = Matrix<int>(n_ + 1, n_, &bsvert_[0][0]);
    cpp_perm_ = Matrix<int>(n_, nsimp_, &perm_[0][0]);

    return;
}

void HyperCube::mkcube(double *cvert_, // cvert_[ncvert_][n_]
        int *bsvert_, // bsvert_[n_ + 1][n_]
        int *perm_, // perm_[n_][nsimp_]
        int ncvert_, int nsimp_, int n_) {
    int i, index, j, k;
    // hypercube vertices
    for (index = 0; index < ncvert_; index++) {
        k = index;
        for (i = 0; i < n_; i++) {
            cvert_[index * n_ + (n_ - 1 - i)] = k - 2 * (k / 2);
            k = k / 2;
        }
    }
    // basic simplex vertices
    for (j = 0; j <= n_; j++) {
        for (i = 0; i < n_; i++) {
            bsvert_[j * n_ + i] = 0;
        }
    }
    for (j = 1; j <= n_; j++) {
        for (i = n_ - j; i < n_; i++) {
            bsvert_[j * n_ + i] = 1;
        }
    }
    // permutations (define  nsimp_)
    if (n_ > 0) {
        mkperm(perm_, n_, nsimp_);
    }
    return;
}

/*
 * This routine creates the array  perm_  containing the permutations
 * of the numbers from 1 to  n_.  It also sets  nsimp_  to be  n_  factorial.
 */
void HyperCube::mkperm(int *perm_, // perm_[n_][nsimp_]
        int n_, int nsimp_) {
    //printf("Inside mkperm. n = %d, nsimp = %d\n", n_, nsimp_);
    int i, j, k, kmfact, nfact, l, shift;
    // initialize
    //    nfact = 0;
    //    if (n_ == 0)   TODO, prestar atencao se no alocamento de espaco eh necesario nsimp = 0!?
    //        return;

    nfact = 1;
    perm_[0] = 0; //TODO, estah certo, pois vai dar os indices de uma matriz!!
    if (n_ <= 1)
        return;
    // loop to build permutations
    for (k = 1; k < n_; k++) {
        kmfact = nfact;
        // place  k  in row  k
        for (j = 0; j < kmfact; j++) {
            perm_[k * nsimp_ + j] = k;
        }
        // place  k  in the other rows
        for (l = 0; l < k; l++) {
            shift = 0;
            for (i = 0; i < k; i++) {
                if (i == k - 1 - l)
                    shift = 1;
                for (j = 0; j < kmfact; j++) {
                    perm_[(i + shift) * nsimp_ + (j + nfact)] = perm_[i * nsimp_ + j];
                }
            }
            for (j = 0; j < kmfact; j++) {
                perm_[(k - 1 - l) * nsimp_ + (j + nfact)] = k;
            }
            nfact = nfact + kmfact;
        }
    }

//    putmi("PERM", perm_, n_, nsimp_);

    return;
}

/**
 * This routine creates the list  face_  of all m_-dimensional faces
 * of the  n_!  simplices obtained by permuting the coordinates of
 * the basic simplex.  to do this it needs the array  comb_  of all
 * subsets of size  m_+1  from a set of size  n_+1.  since  face_  contains
 * only distinct faces, a pointer array  facptr_  is needed to
 * determine which column of  face_  corresponds to each permutation and combination.
 */

int HyperCube::mkface(int *face_, int *facptr_, int *fnbr_, // face_[m_ + 1][dimf_], facptr_[nsimp_][nsface_]
        int dimf_, int nsimp_, int n_, int m_, int nsface_,
        int *bsvert_, int *comb_, // bsvert_[n_ + 1][n_], comb_[numberOfCombinations][m_ + 1]
        int *perm_, int *storn_, int *storm_) {// perm[], storn[n_ + 1], storm[m_ + 1]
    int i, j;
    //printf("Entering: mkfcfp()\n");   // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
    int nface_ = mkfcfp(face_, facptr_, nsface_, dimf_, nsimp_, m_, n_, bsvert_, comb_, perm_, storm_);
    //printf("Aqui nface vale: %d\n", nface_);

    if (nface_ == 0) {
        printf("Error, dimensions do not match. Increase dimf_.\n");
        return ERROR;
    }
    //compute the neighbor array

    //global variable initialization
    //int fnbr_[nsface_][nsface_];

    //zero out  fnbr_
    for (j = 0; j < nsface_; j++) {
        for (i = 0; i < nsface_; i++) {
            fnbr_[j * nsface_ + i] = 0;
        }
    }
    mkfnbr(fnbr_, comb_, n_, m_, nsface_, storn_);

    return nface_;
}

/* Combinatorial     ( n )
                     ( m ),    which means n > m.
   If m < n, it will be done for switching (n, m)  */

int HyperCube::combination(int n, int m) {

    int Up = 1, Bottom = 1;
    if (n > 0 && m > 0) {
        if (n < m) {
            int temp = n;
            n = m;
            m = temp;
        }

        //	    	for(int pont_Up = n; pont_Up > (n-m); pont_Up--) {
        //	    		Up = Up * pont_Up;
        //	    	}
        Up = factorial(n);

        //	    	for(int pont_Bottom = m; pont_Bottom > 1; pont_Bottom--) {
        //	    		Bottom = Bottom * pont_Bottom;
        //	    	}
        Bottom = factorial(m) * factorial(n - m);
    }
    return (Up / Bottom);
}

/** This routine makes the arrays  face_  and  facptr_. It also defines the variable  nface_. */
int HyperCube::mkfcfp(int *face_, int *facptr_, // face_[m_ + 1][dimf_], facptr_[nsimp_][nsface_]
        int nsface_, int dimf_, int nsimp_, int m_, int n_,
        int *bsvert_, int *comb_, // bsvert_[n_ + 1][n_], comb_[numberOfCombinations][m_ + 1]
        int *perm_, int *stor) { // perm_[n_][nsimp_], storm[m_ + 1]


    //local variables
    int c, i, j, index, p, posit, v;

    //loop over all faces to create  face_  and  facptr_
    int nface_ = 0;
    //printf("nsimp_ = %d, nsface_ = %d, m_ = %d, n_ = %d\n", nsimp_, nsface_, m_, n_);
    for (p = 0; p < nsimp_; p++) {
        for (c = 0; c < nsface_; c++) {
            //set indices for each vertex
            for (v = 0; v < m_ + 1; v++) {
                index = 0;
                for (i = 0; i < n_; i++) {
                    //printf("c = %d, m_ = %d, v = %d. c * (m_ + 1) + v = %d.\n\n", c, m_, v, c * (m_ + 1) + v);
                    //printf("n_ = %d, comb_[%d] = %d, n_*comb[%d] = %d.\n\n", n_, c * (m_ + 1) + v, comb_[c * (m_ + 1) + v], c * (m_ + 1) + v, comb_[c * (m_ + 1) + v] * n_);
                    //printf("i = %d, nsimp_ = %d, p = %d, i * nsimp_ + p = %d, perm_[%d] = %d\n\n", i, nsimp_, p, i*nsimp_ + p, i*nsimp_ + p, perm_[i * nsimp_ + p]);
                    //printf("bsvert_[%d] = %d\n\n", comb_[c * (m_ + 1) + v] * n_ + perm_[i * nsimp_ + p], bsvert_[comb_[c * (m_ + 1) + v] * n_ + perm_[i * nsimp_ + p]]);
                    index = 2 * index + bsvert_[comb_[c * (m_ + 1) + v] * n_ + perm_[i * nsimp_ + p]];
                }
                stor[v] = index;
            }
            //store the face_ if it is distinct
            if (nface_ > 0) {
                //printf("Posit antes %d.", posit);
                //posit = search(stor, face_, m_ + 1, dimf_);
                posit = search(stor, face_, m_ + 1, nface_, dimf_);
                //printf("Pase por search, posit = %d, nface = %d\n", posit, nface_);
            } else {
                //printf("Pasamos por posit con %d\n", nface_);
                posit = -1;
            }
            if (posit == -1) {
                nface_ = nface_ + 1;
                if (dimf_ < nface_) {
                    printf("Error the dimension for matrix face is small. [posit = -1].\n");
                    return ERROR;
                }

                //printf("    posit = %d. dimf_ = %d. nface_ = %d\n", posit, dimf_, nface_);

                for (i = 0; i < m_ + 1; i++) {
                    face_[i * dimf_ + (nface_ - 1)] = stor[i];
                }

                posit = nface_ - 1;

                //printf("    After for cycle: posit = %d.\n", posit);
            }
            //set the pointer to the current face_
            //printf("p = %d, nsface_ = %d, c = %d, p * nsface_ + c = %d\n", p, nsface_, c, p * nsface_ + c);
            facptr_[p * nsface_ + c] = posit;
            //printf("After facptr = ...\n");
        }
    }
    //printf("Inside mkfcfp(): nface is %d\n", nface_);
    //printf("dimf_ is %d\n", dimf_);
    return nface_;
}

/**
 * This routine searches the columns of a2  for the column array a1.
 * @param a1 The column to be searched
 * @param a2 The array to be analized
 * @param n1 Number of columns of the array
 * @param n2 Number of lines of the array
 * @returns  <p> -1  ----------> a1  not found in  a2 </p> <p> column no.  --> a1  found in  a2 </p>
 */

int HyperCube::search(int *a1, int *a2, int n1, int n2, int n3) {
    //local variables
    int i, j, search = -1;
    //check for empty  a2
    if (n2 == 0)
        return -1;
    for (j = 0; j < n2; j++) {

        for (i = 0; i < n1; i++) {
            if (a1[i] != a2[i * n3 + j])
                goto lab40;
        }
        search = j;
        return search;
lab40:
        ;
    }
    //fell through.  a1  not found.
    return -1;
}

/*
 * This routine makes the adjacency array  fnbr_  for the faces in
 * the basic simplex. Two faces in comb_ are adjacent if they differ in exactly one vertex.
 */
void HyperCube::mkfnbr(int *fnbr_, int *comb_, //fnbr_[nsface_][nsface_], comb_[numberOfCombinations][m_ + 1]
        int n_, int m_, int nsface_,
        int *stor) { //storn[n_+1]
    //local variables
    int i, j, ndif, nf;




    //loop over the faces
    for (nf = 0; nf < nsface_ - 1; nf++) {
        //set comparison array
        for (i = 0; i < n_ + 1; i++) {
            stor[i] = 0;
        }
        cout << "Valor de n_: " << n_ << " Valor de m_: " << m_ << " Valor de nsface_: " << nsface_ << endl;

        for (i = 0; i < m_ + 1; i++) {
            stor[comb_[nf * (m_ + 1) + i]] = 1;
        }
        //compare with the later faces
        //Panters: na linha que segue acho que estava errado o indice da j e mudamos:
        for (j = nf; j < nsface_; j++) {
            ndif = 0;
            for (i = 0; i < m_ + 1; i++) {
                if (stor[comb_[j * (m_ + 1) + i]] == 0)
                    ndif = ndif + 1;
            }
            if (ndif == 1) {
                fnbr_[nf * nsface_ + j] = 1;
                fnbr_[j * nsface_ + nf] = 1;
            }
        }
    }
    return;
}

/**
 * This routine has two purposes: <ul> <p> (1) Invert the array  facptr_  using two arrays: facept_, ptrf_. </p>
 * <p> (2) Create the attribute arrays facatt_ and ptfatt_. These encode the information about which simplex faces
 * reside on cube faces. </p> </ul>
 */
void HyperCube::mkflst(int *facept_, int *ptrf_, // facept_[faceptInd_][2], ptrf_[nface_][2]
        int *facatt_, int *ptfatt_, // facatt_[fdim_][2], ptfatt_[nface_][2]
        int dimf_, int fdim_, int n_, int m_, int nface_, int nsface_, int nsimp_, int ncvert_, int faceptInd_,
        double *cvert_, int *face_, // cvert_[ncvert_][n_], face_[m_ + 1][dimf_]
        int *facptr_, int *wrk_) { // facptr_[nsimp_][nsface_], wrk_ is storm_[m_ + 1]
    //local variables
    int f, i, next, pos, s, v;
    double sum;
    for (i = 0; i < faceptInd_; i++) {
        facept_[i * 2 + 0] = -1;
        facept_[i * 2 + 1] = -1;
    }
    //invert facptr_
    next = 0;
    for (f = 0; f < nface_; f++) {
        pos = next;
        for (s = 0; s < nsimp_; s++) {
            for (i = 0; i < nsface_; i++) {
                if (facptr_[s * nsface_ + i] == f) {
                    facept_[pos * 2 + 0] = i;
                    facept_[pos * 2 + 1] = s;
                    pos = pos + 1;
                }
            }
        }
        ptrf_[f * 2 + 0] = next;
        ptrf_[f * 2 + 1] = pos - next;
        next = pos;
    }
    //determine face_ attributes
    next = 0;
    for (f = 0; f < nface_; f++) {
        pos = next;
        for (i = 0; i < n_; i++) {
            sum = 0.0;
            for (v = 0; v < m_ + 1; v++) {
                sum = sum + cvert_[(face_[v * dimf_ + f] + 1) * n_ + i];
            }
            if (sum == 0.0) {
                facatt_[pos * 2 + 0] = -(i);
                facatt_[pos * 2 + 1] = mkoppf(f, i, 0, n_, m_, nface_, dimf_, face_, wrk_);
                pos = pos + 1;
            } else if (sum == (double) (m_ + 1)) {
                facatt_[pos * 2 + 0] = +(i);
                facatt_[pos * 2 + 1] = mkoppf(f, i, 1, n_, m_, nface_, dimf_, face_, wrk_);
                pos = pos + 1;
            }
        }
        ptfatt_[f * 2 + 0] = next;
        ptfatt_[f * 2 + 1] = pos - next;
        next = pos;
    }
    return;
}

/**
 * This function determines which face_ is opposite face_  f  in
 * a cube.  It is assumed that face_  f  lies on a cube face_.
 * @param f Face number
 * @param xind Indicates which component is constant.
 * @param value Is that constant value.
 * @param oppvrt Auxiliar vector
 * @returns See the search method return values
 */
int HyperCube::mkoppf(int f, int xind, int value, int n_, int m_, int nface_, int dimf_, int *face_, int *oppvrt) {
    //local variables
    int oppval, sign, v;
    //set the opposite face_ transformation
    if (value == 0)
        sign = +1;
    else
        sign = -1;
    oppval = sign * pow2(n_ - (xind + 1));
    //transform each vertex
    for (v = 0; v < m_ + 1; v++) {
        oppvrt[v] = face_[v * dimf_ + f] + oppval;
    }
    //find the index of the opposite face_
    return search(oppvrt, face_, m_ + 1, nface_, dimf_);
}

/**
 * Funcao para Debug!! Imprime array de inteiros trocando linhas com colunas.
 * <p><h3> Deve ser eliminada apos testes da classe </h3></p>
 * @param name Variable name
 * @param matrix The bidimensional array to be printed
 * @param k1 Size of the array
 * @param k2 Size of the array
 */
/*void printInt(String name, int[] [] matrix, int k1, int k2) {
    int i, j;
    System.out.println(name + ": ");
    for (i = 0; i < k1; i++) {
        for (j = 0; j < k2; j++) {
            System.out.print(matrix[i] [j] + " ");
        }
        System.out.println("");
    }
}
 */
/**
 * Funcao para Debug!! Imprime array de doubles na ordem correta.
 * <p><h3> Deve ser eliminada apos testes da classe </h3></p>
 * @param name Variable name
 * @param matrix The bidimensional array to be printed
 * @param k1 Size of the array
 * @param k2 Size of the array
 */
/*void printDouble(String name, double[] [] matrix, int k1, int k2) {
    int i, j;
    System.out.println(name + ": ");
    for (i = 0; i < k1; i++) {
        for (j = 0; j < k2; j++) {
            System.out.print(matrix[i] [j] + " ");
        }
        System.out.println("");
    }
}
 */
/**
 * Doubles the size of an array
 * @param source Array de origem
 * @return Returns an array containing the elements of source with twice its size
 */
/*Object doubleArray(Object source) {
    int sourceLength = Array.getLength(source);
    Class arrayClass = source.getClass();
    Class componentClass = arrayClass.getComponentType();
    Object result = Array.newInstance(componentClass, sourceLength * 2);
    System.arraycopy(source, 0, result, 0, sourceLength);
    return result;
}
 */


/*AQUI COMECA A PARTE CORRESPONDENTE A CUBSOLVER*/

/*
 * Instituto de Matematica Pura e Aplicada - IMPA
 * Departamento de Dinamica dos Fluidos
 *
 */

// TODO: exstfc eh usado dentro de mksoln, porem ainda nao foi inicializado

int HyperCube::cpp_cubsol(int *solptr_, Matrix<double> &cpp_sol_, int dims_, int *sptr_, int nsoln_,
        double *foncub, int *exstfc, int *face,
        int *facptr_, int dimf_, double *cvert, int ncvert_, int n_, int m_,
        int nsimp_, int nsface_, int nface, double *u, double *g,
        double *x, int *wrki) {

    double sol_[n_][dims_];

    int res = cubsol(solptr_, &sol_[0][0], dims_, sptr_, nsoln_, foncub, exstfc,
            face, facptr_, dimf_, cvert, ncvert_, n_, m_, nsimp_, nsface_, nface, u,
            g, x, wrki);

    cpp_sol_ = Matrix<double>(n_, dims_, &sol_[0][0]);

    return res;
}

int HyperCube::cubsol(int *solptr_, double *sol_, int dims_, int *sptr_, int nsoln_,
        double *foncub, int *exstfc, int *face,
        int *facptr_, int dimf_, double *cvert, int ncvert_, int n_, int m_,
        int nsimp_, int nsface_, int nface, double *u, double *g,
        double *x, int *wrki) {


    // I/O: (int solptr[nsimp_][nsface_], double sol_[n_][dims], int dims_, int sptr[nface_], int nsoln_,
    //       double foncub[m_][ncvert_], int ncvert_, int exstfc[nface], int face[m_ + 1][dimf_],
    //       int facptr_[nsimp_][nsface_], int dimf_, double cvert[ncvert_][n_], int n_, int m_,
    //       int nsimp_, int nsface_, int nface, double u[n_][m_ + 1], double g[m_][m_], double gx[m_],
    //       double x[m_], int wrki[m_])

    // u[n_][m_ + 1] SEEMS to be wrkf1(N,M+1), NOT transposed, as seen in hccube.inc:12.
    // g[m_][m_] could not be mapped to an equivalent in Fortran.
    // x[m_] SEEMS to be wrkf2(M), as seen in hccube.inc:12.
    // wrki[m_] SEEMS to be wrki1(N+1), as seen in hccube.inc:12. Notice that the sizes of the arrays don't match.
    // gx is no longer used. Marchesin & Morante (24-fev-2011).

    //local variables
    //double sol_[n_][dims_];
    //find the solution on each face (define  nsoln_)
//    printf("Entering: mksoln(), nface = %d\n", nface); // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
    nsoln_ = mksoln(sol_, dims_, sptr_, nsoln_, foncub,
            exstfc, face, dimf_, cvert, ncvert_, n_, m_, nface,
            u, g, x, wrki);
    //create the list of solutions for the simplices
    // printf("Entering: smpptr()\n");   // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
    smpptr(solptr_, sptr_, facptr_, nsimp_, nsface_);
    //solptr[nsimp_,nsface_]
    // facptr_[nsimp_][nsface_];
    // sptr_[nface_]????;


    //    cout << "Variaveis inicializadas no contour: " << dims_ << " " << nsoln_ << " " << dimf_ << " " << ncvert_ << " " << n_ << " " << m_ << " " << nsimp_ << " " << nsface_ << " " << nface << endl;
    //    cout << "Valor de nsoln_ : " << nsoln_ <<endl;

    return nsoln_;
}

/**
 * This routine calls  affslv  to obtain the solutions on the faces in  face.
 * @param sol_ Array of solutions
 * @param exstfc Bit array indicating which faces in array "face" are really to be considered by cubsol.
 * Dimension of exstfc should be nface
 */
int HyperCube::mksoln(double *sol_, int dims_, int *sptr_, int nsoln_, double *foncub, int *exstfc, int *face, int dimf_, double *cvert, int ncvert_, int n_, int m_, int nface, double *u, double *g, double *x, int *wrki) {
																																/* //TODO: Ver incongruencia entre sol_[n_][dims] e solution[n_][m_],
																																 aqui o java usaba solution, mas embaixo era preenchido sol_, e parece
																																	que deveria ser sol_[n_][dims], foi o que eu fiz.
																																		    // I/O: (double sol_[n_][dims], int dims_, int sptr[nface_], int nsoln_,
																																	 double foncub[m_][ncvert_], int ncvert_,
																																	    //       int exstfc[nface], int face[m_ + 1][dimf_], int dimf_, 
																																		double cvert[ncvert_][n_], int n_, int m_,
																																		    //       int nface, double u[n_][m_ + 1], double g[m_][m_+1],
																																		 double gx[m_], double x[m_], int wrki[m_]) */
//int contador = 0;
    //local variables
    int flag, i, j, indf, k, v;
//    double g0[m_];
//    double gtemp[m_][m_ + 1];
    double soltemp_[dims_][n_];
    // initializing soltemp_
    for (i = 0; i < dims_; i++) {
        for (j = 0; j < n_; j++) soltemp_[i][j] = -999;
    }
																																	/*
																																	    //printf("dims_ = %d, n_ = %d, nface = %d\n, nsoln_ =%d, dimf_ = %d, m_ =%d,
 																																	nface = %d", dims_, n_, nface, nsoln_, dimf_, m_, nface); */
    // loop over the faces to find each solution (if there is one)
    nsoln_ = -1;
    for (indf = 0; indf < nface; indf++) {
																																		/*       // printf("Inside mksoln(): indf = %d, nface = %d\n", indf, nface); 
																														// ************************ I commented out this line (Morante, Wed 09 Feb 2011 11:17:13 PM BRST). */
        sptr_[indf] = -1;
																															/* //if the indf-face is not to be considered by mksoln, skip it
																											        //        cout <<"Exstfc antes do if "<<exstfc[indf]<<" "<<nsoln_<<endl;
																											        //         printf("Inside mksoln(): exstfc[%d] = %d\n", indf, exstfc[indf]); 
																											// ************************ I commented out this line (Morante, Wed 09 Feb 2011 11:17:13 PM BRST).*/
        if (exstfc[indf] != 0) {
            //set the function values  at the face vertices
            for (k = 0; k < m_ + 1; k++) {
																																	                // nao sei se abaixo eh k ou k+1
                v = face[(k) * dimf_ + indf];
                for (i = 0; i < m_; i++) {
                    g[i * (m_ + 1) + k] = foncub[i * ncvert_ + v];
//                    gtemp[i][k] = g[i * (m_ + 1) + k];
																									                    //printf("v: %d, indf: %d, k: %d, gtemp(i,k): %f\n", v, indf, k, gtemp[i][k]);
                    																												//putmf("FONCUBdentro", foncub, m_, ncvert_);
                }
            }

																														//            putmf("gtemp", &gtemp[0][0], m_, m_ + 1);
																																	//            putmf("gaqui", g, m_, m_ + 1);
            //call the affine solver to solve in the standard simplex
            flag = affslv(x, g, m_, wrki);

																								            // printf("After affslv(): flag = %d\n", flag); // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
																										            //skip the rest if no solution ( note -- pointer initialized to 0 )
																											            //            cout <<"Valor de flag: "<<flag<<endl;
//contador++;
//cout << "In mksoln, flag(" << contador << "/"<< nsoln_+2 << "): " << flag << endl;
//for (int i = 0; i < m_; i++) {
//    for (int j = 0; j < m_+1; j++) {
//        cout << " " << g[i*(m_+1) + j];
//    }
//    cout << endl;
//}
//    cout << endl;

 
            if (flag == 0) {
                //set the pointer to the solution
                nsoln_ = nsoln_ + 1;
                if (dims_ < nsoln_) {
                    printf("Error: Insuficient allocated memory for dims in mksoln routine");
                    return 1;
                    //TODO.... Uma tabela vai ser precisa para isto.... (?)
                    //         Nao ha como saber qual o valor inicial de dims_
                    /*	              dims_ *= 2;
                                        for (j = 0; j < n_; j++) {
                                            sol_[j] = (double[]) doubleArray(sol_[j]);
                                        } */
                }
                // printf("Inside mksol(): indf = %d, nsoln_ = %d\n", indf, nsoln_); // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
                sptr_[indf] = nsoln_;
																										                //                printf("sptr= %d\n", sptr_[indf]);
																											                //                printf("indf= %d\n", indf);
																									                //                printf("nsoln= %d\n", nsoln_);

                //set the face vertices
                for (k = 0; k < m_ + 1; k++) {
                    v = face[k * dimf_ + indf];																				 //printf("v,k,dims %d %d %d\n", v, k, dims_);
                    for (i = 0; i < n_; i++) {
																							                        //TODO: Prestar atencao com a ordem a direita, no cvert!!
                        u[i * (m_ + 1) + k] = cvert[(v) * n_ + i];
																					                        //                        printf("v: %d, indf: %d, k: %d, u(i,k): %f\n", v, indf, k, u[i * (m_ + 1) + k]);
                    }
                }

																														                // soltemp_ eh definida como a trnasposta de sol_
                //transform the solution back from the standard simplex

                afftrn(&soltemp_[nsoln_][0], u, x, n_, m_);
																												                //afftrn (sol_[nsoln_*dims_ + 0], u, x, n_, m_);
            }
        }
    }
    for (i = 0; i < n_; i++) {
        for (j = 0; j < dims_; j++) {
            sol_[i * dims_ + j] = soltemp_[j][i];
        }
    }
    return nsoln_ + 1;
}

//    double auxsol[n_]; NAO PARECE SER USADA!!

/*	    int wrki[m_];
            double u[n_][m_ + 1];
            double g[m_][m_];
            double gx[m_];
            double x[m_];
            double solution[n_][m_];
 */

//ATENCAO em Fortran eram I/O
//    int nface = cFace_.getNumberOfFaces();
//    int[] [] face;
//    face = cFace_.getFaceArray();
//    double[] [] cvert;
//TODO em Fortran era I/O    cvert = cFace_.getVertexCoordinatesArray();
//FECHA ATENCAO

//loop over the faces to find each solution (if there is one)
/*java:    nsoln_ = 0;
    for (indf = 0; indf < nface; indf++) {
        sptr_[indf] = 0;
        //if the indf-face is not to be considered by mksoln, skip it
        if (exstfc[indf] != 0) {
                //boolean zeroFound = false;
            //initialize gx
            v = face[0*dimf_ + indf];
            for (i = 0; i < m_; i++) {
                gx[i] = foncub[i*m_ + v];
            }

            //set the function values  at the face vertices
            for (k = 0; k < m_; k++) {
                v = face[(k + 1)*dimf_ + indf];
                for (i = 0; i < m_; i++) {
                    g[i* m_ +k] = foncub[i*m_ + v + 1]; //TODO: Porque nao falta um +1 apos v?
                }
            }

            //if (zeroFound) {
            //	zeroCount++;
            //}
            //call the affine solver to solve in the standard simplex
            //if(!zeroFound || (zeroCount <= 1)) {
                        flag = affslv(x, gx, g, m_, wrki);
                        //skip the rest if no solution ( note -- pointer initialized to 0 )

                        // solve the system
                        // flag = calculateSolution(x, gx, g, m_, wrki);

                        if (flag == 0) {
                            //set the pointer to the solution
                            nsoln_ = nsoln_ + 1;
                            if (dims_ < nsoln_) {
                                // TODO: Originally in Fortran mksoln returns 1 as an error (?).
                                return 1;
                                //TODO.... Uma tabela vai ser precisa para isto.... (?)
                                //         Nao ha como saber qual o valor inicial de dims_
//	                        dims_ *= 2;
//	                        for (j = 0; j < n_; j++) {
//	                            sol_[j] = (double[]) doubleArray(sol_[j]);

//	                        }
//

// TODO: Parece faltar sptr_[indf] = nsoln_;

                            }
                            //set the face vertices
                            for (k = 0; k < m_ + 1; k++) {
                                v = face[k*dimf_ + indf];
                                for (i = 0; i < n_; i++) {
                                    u[i*(m_+1) + k] = cvert[(v + 1)*n_ + i];//TODO: Prestar atencao com a ordem a direita, no cvert!!
                                }
                            }
                            //transform the solution back from the standard simplex
                            double tmp[n_];
                            tmp[0]=-10;
                            afftrn(&tmp[0], u, x, n_, m_);

                            if(tmp[0] != -10) {
                                sptr_[indf] = nsoln_;
                                int pont_solution;
                                    for(pont_solution = 0; pont_solution < n_; pont_solution++) {
                                        sol_[pont_solution*dims_ + (nsoln_ - 1)] = tmp[pont_solution];
                                    }
                            } else {
                                // TODO: In the same vein, we return 1 as error.
                                return 1;
                            }
                        }
            //}

        }
    }
    return 0;
}java:*/

/**
 * This routine finds the zero  x  of the affine function whose values at
 * the vertices of the standard m_-simplex are stored in the columns of  g.
 * It also checks whether the solution is inside the simplex and sets flag appropriately.
 * @param x
 * @param g0
 * @param g
 * @param wrki
 * @return  <ul> <p> -1      Error in solving the system </p> <p>  0      Valid solution inside simplex </p>
 * <p>  1      Valid solution outside simplex </p> <ul>
 */

/* void affslv(double[] x, double[] g0, double[] [] g, int[] wrki) {
 //local variables
 int i, j;
 //subtract  g0  from  g1, g2, ..., gm
 for (j = 0; j < m_; j++) {
     for (i = 0; i < m_; i++) {
         g[i] [j] = g[i] [j] - g0[i];
     }
 }
 //change sign of  g0
 for (i = 0; i < m_; i++) {
     g0[i] = -g0[i];
 }

}*/

int HyperCube::affslv(double *x, double *g, int m_, int *wrki) {
    // x[m_], g[m_][m_+1], m_, wrki[m_]
    //local variables
    double sum;
    int i, j, result;
    double gt[m_][m_];
    double g0[m_];
    //subtract  g0  from  g1, g2, ..., gm
    for (j = 1; j < m_ + 1; j++) {
        for (i = 0; i < m_; i++) {
            g[i * (m_ + 1) + j] = g[i * (m_ + 1) + j] - g[i * (m_ + 1) + 0];
            // TODO: Ver se g esta ao direito ou o inves. (transposta??)
        }
    }
    //change sign of  g0
    for (i = 0; i < m_; i++) {
        g[i * (m_ + 1) + 0] = -g[i * (m_ + 1) + 0];
    }

    //temporal vectors gt and g0

    for (i = 0; i < m_; i++) {
        for (j = 0; j < m_; j++) {
            gt[i][j] = g[i * (m_ + 1) + j + 1];
            // TODO: Ver se g esta ao direito ou o inves. (transposta??)
        }
        g0[i] = g[i * (m_ + 1) + 0];
    }
    //    putmf("ginside", g, m_, m_+1);
    //solve the system
    result = solver(m_, &gt[0][0], g0, x); //TODO: De incluir ipiv[n] aqui deve entrar ", wrki);"
    //check for no solution

//cout << " x:" << x[0] << " " << x[1] << " " << x[2] << endl;

    if (result != 0)
        return result;
    //determine whether solution is inside
    sum = 0.0;
    for (i = 0; i < m_; i++) {
    /* This NEW condition is made in order to exclude repeated points at common edges */
//        if ( (x[i] < 0.0 + KLUDGE) || (x[i] > 1.0 - KLUDGE) ) {
        if (x[i] <= 0.0) {
            //solution is outside
            return 1;
        } else {
            sum = sum + x[i];
        }
    }
    /* This NEW condition is made in order to exclude repeated points at common edges */
//    if ( (sum <= 0.0 + KLUDGE) || (sum >= 1.0 - KLUDGE) ) {
    if (sum >= 1.0) {
        //solution is outside
        return 1;
    }
    //fell through -- solution is inside
    return 0;
}

/**
 * <p> This routine maps  x  in the standard m_-simplex to  y  in the
 * n_-dimensional m_-simplex with vertices  u(.,1), u(.,2), ...,  u(.,m_+1). </p>
 * <p> The map is the affine transformation which takes the vertices of the
 * standard simplex to those contained in the columns of  u. </p>
 * @param y
 * @param u
 * @param x
 */
void HyperCube::afftrn(double *solution, double *u, double *x, int n_, int m_) {
    // solution[n_], u[n_][m_ + 1], x[m_], n_, m_
    int i, k;
    //    printf("Entrando en afftrn\n");
    for (i = 0; i < n_; i++) {
        //        printf("solution[%d] = %f\n", i, solution[i]);
        //        printf("u[%d*(m_ + 1) + 0] = %f\n", i, u[i * (m_ + 1) + 0]);
        solution[i] = u[i * (m_ + 1) + 0];
        for (k = 0; k < m_; k++) {
            solution[i] += (u[i * (m_ + 1) + (k + 1)] - u[i * (m_ + 1) + 0]) * x[k];
            /*if (sol_[i] [ind] == Double.NEGATIVE_INFINITY) {
                sol_[i] [ind] = -Double.MIN_VALUE;
            } else if (sol_[i] [ind] == Double.POSITIVE_INFINITY) {
                sol_[i] [ind] = Double.MAX_VALUE;
            } */
        }
    }
    return;
}

int HyperCube::solver(int n, double *A, double *b, double *x) {
    /**
     * Seems that ipiv[n] is missing, but with the new solver it is not necessary to keep the
     * line permutation (index pivoting).
     **/

    double eps = 1e-12;

    double det, anorm;
    int i, j;

    if (n == 1) {
        if (fabs(A[0]) <= eps) return -1;
        x[0] = b[0] / A[0];

        return 0;
    } else if (n == 2) {
        det = A[0] * A[3] - A[1] * A[2];
        anorm = 0;
        for (i = 0; i < n * n; i++) anorm += A[i] * A[i];

        if (fabs(det) <= (eps * anorm)) return -1;

        x[0] = (b[0] * A[3] - b[1] * A[1]) / det;
        x[1] = (A[0] * b[1] - A[2] * b[0]) / det;

        return 0;
    } else if (n == 3) {
        det = A[0]*(A[4] * A[8] - A[5] * A[7])
                - A[3]*(A[1] * A[8] - A[7] * A[2])
                + A[6]*(A[1] * A[5] - A[4] * A[2]);

        anorm = 0;
        for (i = 0; i < n * n; i++) anorm += A[i] * A[i];

        if (fabs(det) <= (eps * anorm)) return -1;

        x[0] = (b[0]*(A[4] * A[8] - A[7] * A[5])
                - b[1]*(A[1] * A[8] - A[7] * A[2])
                + b[2]*(A[1] * A[5] - A[4] * A[2])
                ) / det;

        x[1] = (-b[0] * (A[3] * A[8] - A[6] * A[5])
                + b[1] * (A[0] * A[8] - A[6] * A[2])
                - b[2] * (A[0] * A[5] - A[3] * A[2])
                ) / det;

        x[2] = (b[0] * (A[3] * A[7] - A[6] * A[4])
                - b[1] * (A[0] * A[7] - A[6] * A[1])
                + b[2] * (A[0] * A[4] - A[3] * A[1])
                ) / det;
        return 0;
    } else {
        int dim = n;
        int nrhs = 1;

        int lda = n;
        int ipiv[n];
        int ldb = n;
        int info;

        // Create a transposed copy of A to be used by LAPACK's dgesv:
        double B[n][n];
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) B[j][i] = A[i * n + j];
        }

        // Create a copy of b to be used by LAPACK's dgesv:
        double bb[n];
        for (i = 0; i < n; i++) bb[i] = b[i];

        dgesv_(&dim, &nrhs, &B[0][0], &lda, &ipiv[0], &bb[0], &ldb, &info);

        if (info == 0) {
            for (i = 0; i < n; i++) x[i] = bb[i];
            return 0;
        } else return -1;
    }
}

void HyperCube::smpptr(int *solptr_, int *sptr_, int *facptr, int nsimp_, int nsface_) {
    //solptr[nsimp_,nsface_]
    // facptr_[nsimp_][nsface_];
    //sptr_[nface_]????
    //local variables
    int ns, nf;
    //loop over  facptr  to find solutions
    for (ns = 0; ns < nsimp_; ns++) {
        for (nf = 0; nf < nsface_; nf++) {
            solptr_[ns * nsface_ + nf] = sptr_[facptr[ns * nsface_ + nf]];
            // printf("Inside smpptr(): nsface = %d, solptr[%d][%d] = %d\n", nsface_, ns, nf, solptr_[nf]); // I commented out this line (Morante: Wed 09 Feb 2011 11:20:26 PM BRST).
        }
    }
    return;
}


/** The following routine builds the array of pointers to the edges_ in the solution polyhedron. */

/*input
   -----
     solptr  --  a 2-dimensional array containing pointers to the
                 solution in  sol.  i.e.,  solptr(i,j)  is the
                 column of  sol  containing the solution for the
                 i-th face in the j-th simplex.  (if  solptr  is
                 0, the corresponding face contains no solution.)

     fnbr    --  a 2-dimensional array containing the adjacency
                 information for the faces in the basic simplex.
                 i.e.,  fnbr(i,j)  is 1 (resp. 0) if the i-th and
                 j-th faces are (resp. are not) adjacent in the
                 basic simplex.  here,  i,j  refer to the order
                 determined by the array  comb.

     nsimp   --  the number of simplices in the hypercube.

     nsface  --  the number of faces in the basic simplex.

   output
   ------
     edges   --  a 2-dimensional array each column of which contains
                 pointers to the array  sol  at which the vertices
                 of an edge of the solution polyhedron are stored.

     nedges  --  the number of edges found.

     smpedg  --  a 2-dimensional array each column of which contains
                 the beginning and ending indices for the array
                 edges.  i.e., the edges for simplex number  ns
                 are in columns  smpedg(1,ns)  to  smpedg(2,ns)
                 of  edges.  there are no edges if
                        smpedg(1,ns) .gt. smpedg(2,ns).
c*/

int HyperCube::cpp_mkedge(Matrix<int> &cpp_edges_, int dime_, int nedges_, int *smpedg_, int *solptr_, int *fnbr_, int nsimp_, int nsface_) {
    int edges_[2][dime_];

    int res = mkedge(&edges_[0][0], dime_, nedges_, smpedg_, solptr_, fnbr_, nsimp_, nsface_);

    cpp_edges_ = Matrix<int>(2, dime_, &edges_[0][0]);

    return res;
}

int HyperCube::mkedge(int *edges_, int dime_, int nedges_, int *smpedg_, int *solptr_, int *fnbr_, int nsimp_, int nsface_) {

    //solptr[nsimp_,nsface_]
    //fnbr_[nsface_][nsface_]
    //edges_[2][dime_], smpedg_[nsimp_][2]

    //Global variables initialization
    //int edges_[2] [dime_];
    // int smpedg_[nsimp] [2];
    //local variables
    int i, j, ns, spi, spj;
    //initialize
    nedges_ = -1;
    //loop over the simplices creating the edges_
    for (ns = 0; ns < nsimp_; ns++) {
        smpedg_[ns * 2 + 0] = nedges_ + 1;
        //determine which neighboring faces have a solution edge
        for (i = 0; i < nsface_ - 1; i++) {
            // printf("ns = %d, nsface_ = %d, i = %d\n", ns, nsface_, i);  // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
            spi = solptr_[ns * nsface_ + i]; // printf("spi = %d\n", spi);  // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
            if (spi != -1) {
                for (j = i + 1; j < nsface_; j++) {
                    if (fnbr_[j * nsface_ + i] != 0) {
                        // printf("\nns = %d, nsface_ = %d, j = %d\n", ns, nsface_, j);  // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
                        spj = solptr_[ns * nsface_ + j]; // printf("spj = %d\n", spj);   // I commented out this line (Morante: Wed 09 Feb 2011 11:17:50 PM BRST )
                        if (spj != -1) {
                            if (dime_ < nedges_ + 1) {
                                printf("Error: Insuficient allocated memory for dime in mkedge routine");
                                smpedg_[ns * 2 + 1] = nedges_;
                                return 1;
                            }
                            nedges_ = nedges_ + 1;
                            edges_[0 * dime_ + nedges_] = spi;
                            edges_[1 * dime_ + nedges_] = spj;
                        }
                    }
                }
            }
        }
        smpedg_[ns * 2 + 1] = nedges_;
    }
    return nedges_ + 1;
}

/**
 *  This routine makes a list of all  m_  element subsets of an n_  element set.  It also sets
 *  nsface_  to be the number of subsets found. Each subset is stored in a column of  comb_.
 *  The elements in a column are in increasing order, and the columns are ordered lexicographically.
 *  @param n_ Set size
 *  @param m_ Subset size
 */
int HyperCube::mkcomb(int *comb_, // comb_[numberOfCombinations][m_ + 1]
        int np, int mp) { // np = n_ + 1; mp = m_ + 1
    // notice that m is equal to m_+1
    int i, j, pos, count;
    //    int aux = m;

    int nsface_ = 1;
    for (i = 0; i < mp; i++) {
        comb_[(nsface_ - 1) * mp + i] = i;
    }
    //the main loop
    while (1) {
        count = 0;
        //search for the last position to change
        for (i = 0; i < mp; i++) {
            pos = mp - 1 - i;
            if (comb_[(nsface_ - 1) * mp + pos] == np - i - 1) {
                count++;
                //return;
            } else {
                //position found.  store new combination.
                nsface_++;
                //increment the value at the position
                comb_[(nsface_ - 1) * mp + pos] = comb_[(nsface_ - 2) * mp + pos] + 1;
                //copy the preceding values
                for (j = 0; j < pos; j++) {
                    comb_[(nsface_ - 1) * mp + j] = comb_[(nsface_ - 2) * mp + j];
                }
                //modify the succeeding values
                for (j = pos + 1; j < mp; j++) {
                    comb_[(nsface_ - 1) * mp + j] = comb_[(nsface_ - 1) * mp + (j - 1)] + 1;
                }
                //break
                i = mp;
            }
        }
        if (count == mp) {
            //fell through.  done.
            return nsface_;
        }
    }
    return nsface_;
}

/*
 * Find facept_ and facatt_ sizes
 */
void HyperCube::findInds(int *faceptInd_, int *fdim_, //pointers to initializate these values
        int n_, int m_,
        double *cvert_, int *face_, // cvert_[ncvert_][n_], face_[m_ + 1][dimf_]
        int dimf_, int nface_, int nsimp_, int nsface_, int *facptr_) {
    //local variables
    int f, i, pos, s, v;
    double sum;
    pos = 0;
    for (f = 0; f < nface_; f++) {
        for (s = 0; s < nsimp_; s++) {
            for (i = 0; i < nsface_; i++) {
                if (facptr_[s * nsface_ + i] == f) {
                    pos = pos + 1;
                }
            }
        }
    }
    //Global variable initialization
    *faceptInd_ = pos + 1;
    pos = 0;
    for (f = 0; f < nface_; f++) {
        for (i = 0; i < n_; i++) {
            sum = 0.0;
            for (v = 0; v < m_ + 1; v++) {
                sum = sum + cvert_[(face_[v * dimf_ + f]) * n_ + i];
            }
            if (sum == 0.0) {
                pos = pos + 1;
            } else if (sum == (double) (m_ + 1)) {
                pos = pos + 1;
            }
        }
    }
    *fdim_ = pos;
    return;
}

/**
     the following routine builds the array of pointers to the edges in
c     the solution polyhedron which are parallel to a coordinate plane.
c
c   input
c   -----
c     solptr  --  a 2-dimensional array containing pointers to the
c                 solution in  sol.  i.e.,  solptr(i,j)  is the
c                 column of  sol  containing the solution for the
c                 i-th face in the j-th simplex.  (if  solptr  is
c                 0, the corresponding face contains no solution.)
c
c     fnbr    --  a 2-dimensional array containing the adjacency
c                 information for the faces in the basic simplex.
c                 i.e.,  fnbr(i,j)  is 1 (resp. 0) if the i-th and
c                 j-th faces are (resp. are not) adjacent in the
c                 basic simplex.  here,  i,j  refer to the order
c                 determined by the array  comb.
c
c     nsimp   --  the number of simplices in the hypercube.
c
c     nsface  --  the number of faces in the basic simplex.
c
c   output
c   ------
c     edges   --  a 2-dimensional array each column of which contains
c                 pointers to the array  sol  at which the vertices
c                 of an edge of the solution polyhedron are stored.
c
c     nedges  --  the number of edges found.
c
c     smpedg  --  a 2-dimensional array each column of which contains
c                 the beginning and ending indices for the array
c                 edges.  i.e., the edges for simplex number  ns
c                 are in columns  smpedg(1,ns)  to  smpedg(2,ns)
c                 of  edges.  there are no edges if
c                        smpedg(1,ns) .gt. smpedg(2,ns).
c
 */
int HyperCube::mklevl(int *edges_, int dime_, int nedges_, int *smpedg_, int *solptr_, int *fnbr_, int nsimp, int nsface, int *facatt, int *ptfatt, int *facptr, int fdim_, int nface_) {

    // edges_[2][dime_], smpedg_[nsimp_][2], solptr[nsimp_,nsface_]
    //fnbr_[nsface_][nsface_], facptr_[nsimp_][nsface_]
    // facatt_[fdim_][2], ptfatt_[nface_][2]
    //edges_[2][dime_], smpedg_[nsimp_][2]

    //    putmi("facatt", facatt, fdim_, 2);
    //    putmi("ptfatt", ptfatt, nface_, 2);



    //local variables
    int i, j, ns, spi, spj;
    int begi, begj, k, l;

    //initialize
    nedges_ = 0;
    //loop over the simplices creating the edges_
    for (ns = 0; ns < nsimp; ns++) {
        smpedg_[ns * 2 + 0] = nedges_ + 1;
        //determine which neighboring faces have a solution edge
        for (i = 0; i < nsface - 1; i++) {
            spi = solptr_[ns * nsface + i];
            if (spi != -1) {
                for (j = i + 1; j < nsface; j++) {
                    if (fnbr_[j * nsface + i] != 0) {
                        spj = solptr_[ns * nsface + j];
                        if (spj != -1) {
                            begi = ptfatt[(facptr[ns * nsface + i])*2 + 0];
                            begj = ptfatt[(facptr[ns * nsface + j])*2 + 0];
                            for (k = 0; k < ptfatt[(facptr[ns * nsface + i])*2 + 1]; k++) {
                                for (l = 0; l < ptfatt[(facptr[ns * nsface + j]) *2 + 1]; l++) {
                                    if (facatt[(begi + k - 1) * 2 + 0] != facatt[(begj + l - 1) * 2 + 0]) continue;
                                    else goto label40;
                                }
                                goto label60;
label40:
                                ;
                            }

                            //found an edge
                            if (dime_ < nedges_ + 1) {
                                printf("Error: Insuficient allocated memory for dime in mklevl routine");
                                smpedg_[ns * 2 + 1 ] = nedges_;
                                return 1;
                            }
                            nedges_ = nedges_ + 1;
                            //printf("PASAMOS POR AQUI CON spi, spj = %d, %d", spi, spj);
                            edges_[0 * dime_ + nedges_] = spi;
                            edges_[1 * dime_ + nedges_] = spj;
                        }
                    }
label60:
                    ;
                }
            }
        }

        smpedg_[ns * 2 + 1] = nedges_;
    }
    return nedges_ + 1;
}





// The following routine prints a real matrix, priorizing the longest lenght as lines.
// For traking the actual size is printed within the name

void HyperCube::putmf(const char *name, double *matrix_, int k1_, int k2_) {
    int i, j, auxtemp = 0;
    int IndFirst = 0, IndLast;
    printf("DOUBLE ARRAY -- \t %s[%d, %d]:\n", name, k1_, k2_);

    if (k1_ > k2_) {
        auxtemp = k1_;
        k1_ = k2_;
        k2_ = auxtemp;
    }

    IndLast = k2_;

    do {
        if (IndLast - IndFirst > 8) IndLast = IndFirst + 8;

        for (j = IndFirst; j < IndLast; j++) printf("%9d \t", j);
        printf("\n");

        for (j = IndFirst; j < IndLast; j++) printf("----------\t");
        printf("\n");

        if (auxtemp == 0) {
            for (i = 0; i < k1_; i++) {
                for (j = IndFirst; j < IndLast - 1; j++) {
                    printf("% 3.2E    \t", matrix_[i * k2_ + j]);
                }
                printf("% 3.2E\n", matrix_[i * k2_ + IndLast - 1]);
            }
        } else {
            for (i = 0; i < k1_; i++) {
                for (j = IndFirst; j < IndLast - 1; j++) {
                    printf("% 3.2E    \t", matrix_[i + j * k1_]);
                }
                printf("% 3.2E\n", matrix_[i + (IndLast - 1) * k1_]);
            }
        }
        printf("\n");

        IndFirst = IndLast;
        IndLast = k2_;

    } while (IndLast - IndFirst > 0);

    return;
}

void HyperCube::putmf2(const char *name, double *matrix_, int k1_, int k2_, int kmax) {
    int i, j;
    double temporal_matrix[k1_][kmax];

    for (i = 0; i < k1_; i++) {
        for (j = 0; j < kmax; j++) {
            temporal_matrix[i][j] = matrix_[i * k2_ + j];
        }
    }

    putmf(name, &temporal_matrix[0][0], k1_, kmax);

    return;
}




// The following routine prints an integer matrix, priorizing the longest lenght as lines.
// For traking the actual size is printed within the name

void HyperCube::putmi(const char *name, int *matrix_, int k1_, int k2_) {
    int i, j, auxtemp = 0;
    int IndFirst = 0, IndLast;
    printf("INTEGER ARRAY -- \t %s[%d, %d]:\n", name, k1_, k2_);

    if (k1_ > k2_) {
        auxtemp = k1_;
        k1_ = k2_;
        k2_ = auxtemp;
    }

    IndLast = k2_;

    do {
        if (IndLast - IndFirst > 30) IndLast = IndFirst + 30;

        for (j = IndFirst; j < IndLast; j++) printf("%3d ", j);
        printf("\n");

        for (j = IndFirst; j < IndLast; j++) printf("--- ");
        printf("\n");

        if (auxtemp == 0) {
            for (i = 0; i < k1_; i++) {
                for (j = IndFirst; j < IndLast; j++) {
                    printf("%3d ", matrix_[i * k2_ + j]);
                }
                printf("\n");
            }
        } else {
            for (i = 0; i < k1_; i++) {
                for (j = IndFirst; j < IndLast; j++) {
                    printf("%3d ", matrix_[i + j * k1_]);
                }
                printf("\n");
            }
        }
        printf("\n");

        IndFirst = IndLast;
        IndLast = k2_;

    } while (IndLast - IndFirst > 0);

    return;
}

void HyperCube::putmi2(const char *name, int *matrix_, int k1_, int k2_, int kmax) {
    int i, j;
    int temporal_matrix[k1_][kmax];

    for (i = 0; i < k1_; i++) {
        for (j = 0; j < kmax; j++) {
            temporal_matrix[i][j] = matrix_[i * k2_ + j];
        }
    }

    putmi(name, &temporal_matrix[0][0], k1_, kmax);

    return;
}

//int mkpoly( int *polyv_, int dimv, int npolys, int *smppol_, int *solptr_, int *fnbr_, int nsimp_, int nsface_, int dimp_) {
//
//    //  smppol_[dimp_][2], solptr[nsimp_,nsface_]
//            //fnbr_[nsface_][nsface_], f
//            //polyv_[dimv_],
//
//    // local variables
//
//    int cur, i, l, ns, nverts, nvertt, prev, solcnt, spi;
//
//    //initialize
//    npolys=0;
//    nverts=0;
//
//    //loop over the simplices creating the polygons
//
//    for (ns=0; ns < nsimp_; ns++){
//    //set a temporary vertex count
//        nvertt=nverts;
//        //find the first solution in this simplex and count the number
//        solcnt=0;
//        for (i=0; i < nsface_; i++){
//            spi=solptr_[ns*nsface_ +i];
//            if (spi != 0){
//                if (solcnt == 0){
//                    if (dimv <= nvertt){
//                        printf("Error in mkpoly: insufficient dimv");
//                        return 2;
//                    }
//                    nvertt++;
//                    polyv_[nvertt] = spi;
//                    prev = 0;
//                    cur = i;
//                }
//                solcnt=solcnt+1;
//            }
//        }
//        // check for no polygon
//        if (solcnt > 0){
//            for (l=1; l < solcnt; l++){
//                for (i=0; i < nsface_; i++){
//                    if (i == prev || i==cur) continue;
//                    spi = solptr_[ns*nsface_+i];
//                    if (spi != 0 && fnbr_[i*nsface_ + cur] != 0){
//                        if (dimv <= nvertt){
//                            printf("Error in mkpoly: increase dimv");
//                            return 2;
//                        }
//                        nvertt++;
//                        polyv_[nvertt] = spi;
//                        prev = cur;
//                        cur = i;
//                        goto label140;
//                    }
//                }
//
//                // error:polygon closed too soon
//                printf("Error in mkpoly: partial polygon");
//                goto label200;
//                label140:;
//            }
//
//            if (dimp_ <= npolys){
//                printf("Error in mkpoly: insufficient dimp");
//                return 1;
//            }
//            npolys++;
//            smppol_[npolys *2 + 0] = nvertt + 1;
//            smppol_[npolys *2 + 1] = nvertt;
//            nverts = nvertt;
//        }
//        label200:;
//    }
//    return 0;
//}

// The following routine prints an integer matrix, priorizing the longest lenght as lines.
// For traking the actual size is printed within the name

//void putmi(const char *name, int *matrix_, int k1_, int k2_) {
//    int i, j, auxtemp = 0;
//    int IndFirst = 0, IndLast;
//    printf("INTEGER ARRAY -- \t %s[%d, %d]:\n", name, k1_, k2_);
//
//    if (k1_ > k2_) {
//	auxtemp = k1_;
//        k1_ = k2_;
//        k2_ = auxtemp;
//    }
//
//    IndLast = k2_;
//
//    do {
//        if (IndLast - IndFirst > 30) IndLast = IndFirst + 30;
//
//        for (j = IndFirst; j < IndLast; j++) printf("%3d ", j);
//        printf("\n");
//
//        for (j = IndFirst; j < IndLast; j++) printf("--- ");
//        printf("\n");
//
//        if (auxtemp == 0) {
//            for (i = 0; i < k1_; i++) {
//                for (j = IndFirst; j < IndLast; j++) {
//                    printf("%3d ", matrix_[i*k2_ + j]);
//                }
//                printf("\n");
//            }
//        } else {
//            for (i = 0; i < k1_; i++) {
//                for (j = IndFirst; j < IndLast; j++) {
//                    printf("%3d ", matrix_[i + j*k1_]);
//                }
//                printf("\n");
//            }
//        }
//        printf("\n");
//
//	IndFirst = IndLast;
//        IndLast = k2_;
//
//    } while (IndLast - IndFirst > 0);
//
//    return;
//}

//void putmi2(const char *name, int *matrix_, int k1_, int k2_, int kmax) {
//    int i, j;
//    int temporal_matrix[k1_][kmax];
//
//    for (i = 0; i < k1_; i++) {
//	for (j = 0; j < kmax; j++) {
//	    temporal_matrix[i][j] = matrix_[i*k2_ + j];
//	}
//    }
//
//    putmi(name, &temporal_matrix[0][0], k1_, kmax);
//
//    return;
//}

// The following routine prints a real matrix, priorizing the longest lenght as lines.
// For traking the actual size is printed within the name

//void putmf(const char *name, double *matrix_, int k1_, int k2_) {
//    int i, j, auxtemp = 0;
//    int IndFirst = 0, IndLast;
//    printf("DOUBLE ARRAY -- \t %s[%d, %d]:\n", name, k1_, k2_);
//
//    if (k1_ > k2_) {
//	auxtemp = k1_;
//        k1_ = k2_;
//        k2_ = auxtemp;
//    }
//
//    IndLast = k2_;
//
//    do {
//        if (IndLast - IndFirst > 8) IndLast = IndFirst + 8;
//
//        for (j = IndFirst; j < IndLast; j++) printf("%9d \t", j);
//        printf("\n");
//
//        for (j = IndFirst; j < IndLast; j++) printf("----------\t");
//        printf("\n");
//
//        if (auxtemp == 0) {
//            for (i = 0; i < k1_; i++) {
//                for (j = IndFirst; j < IndLast - 1; j++) {
//                    printf("% 3.2E    \t", matrix_[i*k2_ + j]);
//                }
//                printf("% 3.2E\n", matrix_[i*k2_ + IndLast - 1]);
//            }
//        } else {
//            for (i = 0; i < k1_; i++) {
//                for (j = IndFirst; j < IndLast - 1; j++) {
//                    printf("% 3.2E    \t", matrix_[i + j*k1_]);
//                }
//                printf("% 3.2E\n", matrix_[i + (IndLast - 1)*k1_]);
//            }
//        }
//        printf("\n");
//
//	IndFirst = IndLast;
//        IndLast = k2_;
//
//    } while (IndLast - IndFirst > 0);
//
//    return;
//}

//void putmf2(const char *name, double *matrix_, int k1_, int k2_, int kmax) {
//    int i, j;
//    double temporal_matrix[k1_][kmax];
//
//    for (i = 0; i < k1_; i++) {
//	for (j = 0; j < kmax; j++) {
//	    temporal_matrix[i][j] = matrix_[i*k2_ + j];
//	}
//    }
//
//    putmf(name, &temporal_matrix[0][0], k1_, kmax);
//
//    return;
//}

//void putme(const char *name, double *matrix_, int k1_, int k2_) {
//    int i, j, k;
//    printf("DOUBLE ARRAY -- \t %s[%d, %d, 2]:\n", name, k1_, k2_);
//
//    for (j = 0; j < k1_; j++) printf("%9d ", j);
//    printf("\n");
//
//    for (j = 0; j < k1_; j++) printf("----------\t");
//    printf("\n");
//
//    for (j = 0; j < k2_; j++) {
//          for (i = 0; i < k1_; i++) {
//            for (k =0; k < 2; k++){
//           printf("% 3.2E    \t", matrix_[(2 * i + k)*k2_ + j]);
//           }
//           printf("\n");
//        }
//    }
//    printf("\n");
//    return;
//}

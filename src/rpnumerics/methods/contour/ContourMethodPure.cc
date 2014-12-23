/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) HugoniotContourMethodPure.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "ContourMethodPure.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

// Static variables are defined here:
//
bool ContourMethodPure::is_first = true;

HyperCube ContourMethodPure::hc;

int ContourMethodPure::hn = 2;
int ContourMethodPure::hm = 1;

int ContourMethodPure::nsface_ = 3;
int ContourMethodPure::nface_ = 5;
int ContourMethodPure::nsoln_ = -1;
int ContourMethodPure::nedges_;

int ContourMethodPure::dims_ = 6; // Estava 50, Fortran diz 6 (inicialmente 4)
int ContourMethodPure::dime_ = 6; // Estava 60, Fortran diz 6 (inicialmente 2)
int ContourMethodPure::dimf_ = 5;
int ContourMethodPure::ncvert_ = 4;
int ContourMethodPure::nsimp_ = 2;

int ContourMethodPure::numberOfCombinations = 3;

int * ContourMethodPure::storn_;
int * ContourMethodPure::storm_;

double * ContourMethodPure::cvert_;
double * ContourMethodPure::vert;
int * ContourMethodPure::bsvert_;
int * ContourMethodPure::perm_;
int * ContourMethodPure::comb_;
int * ContourMethodPure::fnbr_;
int * ContourMethodPure::facptr_;
int * ContourMethodPure::face_;
double * ContourMethodPure::sol_;
int * ContourMethodPure::solptr_;
int * ContourMethodPure::edges_;
int * ContourMethodPure::smpedg_;
int * ContourMethodPure::exstfc;
int * ContourMethodPure::sptr_;

int ContourMethodPure::tsimp = 1;
int ContourMethodPure::tface = 3;

void ContourMethodPure::allocate_arrays(void){
    if (is_first){
        storn_ = new int[hn + 1];
        storm_ = new int[hm + 1];

        cvert_ = new double[ncvert_*hn];
        vert = new double[ncvert_*hn];
        bsvert_ = new int[(hn + 1)*hn];
        perm_ = new int[hn*nsimp_];
        comb_ = new int[numberOfCombinations*(hm + 1)];

        nsface_ = hc.mkcomb(comb_, hn + 1, hm + 1);
        fnbr_ = new int[nsface_*nsface_];

        //inicializing another arrays, it were globally defined in java
        facptr_ = new int[nsimp_*nsface_];
        face_ = new int[(hm + 1)*dimf_];

        hc.mkcube(cvert_, bsvert_, perm_, ncvert_, nsimp_, hn);
        nface_ = hc.mkface(face_, facptr_, fnbr_, dimf_, nsimp_, hn, hm, nsface_,
                           bsvert_, comb_, perm_, storn_, storm_);

        sol_ = new double[hn*dims_];
        solptr_ = new int[nsimp_*nsface_];
        edges_ = new int[2*dime_];
        smpedg_ = new int[nsimp_*2];

        exstfc = new int[nface_];
        sptr_ = new int[nface_];

        is_first = false;

        printf("++++++++++++++++ REMEMBER TO INVOKE deallocate_arrays() AT QUIT-TIME!!!\n++++++++++++++++ DON\'T SAY I DIDN\'T WARN YOU!!!\n");
    }

    return;
}

void ContourMethodPure::deallocate_arrays(void){
    if (!is_first){
        delete sptr_;
        delete exstfc;

        delete smpedg_;
        delete edges_;
        delete solptr_;
        delete sol_;

        delete face_;
        delete facptr_;

        delete fnbr_;

        delete comb_;
        delete perm_;
        delete bsvert_;
        delete vert;
        delete cvert_;

        delete storm_;
        delete storn_;

        is_first = true;
    }

    return;
}

ContourMethodPure::ContourMethodPure(HugoniotFunctionClass *h){
    hugoniot = h;
}

ContourMethodPure::~ContourMethodPure(){
}

//double ContourMethodPure::f(double x, double y) {
////    return 3.0*x-y-1.2;//pow(x - 0.3, 2.0) + pow(y - 0.7, 2.0) - 0.01;

//    // return pow(x - 0.313, 2.0) + pow(y - 0.724, 2.0) - 0.017;
//    return pow(x*x + y*y, 2.0) - 2.0*(x*x - y*y);
//    //    return 1 * (x - 0.1) + 2 * (y - 0.2);
//}

//int ContourMethodPure::inpdom(double *u) { // double u[2]//Replace by Boundary::inside
////        if (u[0] >= 0 && u[1] >= 0 && u[0] + u[1] <= 1) return 1;
////        else return 0;

//     return 1;
//}

//int ContourMethodPure::curv2d(/*double *segend,*/ int sn, int seglim, double fdummy, double *rect, int *res, int ifirst) {
int ContourMethodPure::curv2d(/*double *segend,*/ int sn, int seglim, double fdummy, double *rect, int *res, int ifirst,
                          std::vector<RealVector> &vrs) { //TODO: int sn must be eliminated (because it is == vrs.size()).
                                                          // TODO: Get rid of seglim & fdummy, ifirst (use ctor instead), and rect will become Domain*.
//    printf("BEGINS: curv2d()\n");
////    HyperCube hc;
//    vrs.clear();

////    double segend[seglim][2][2]; //int sn, int seglim, double f, double rect[4], int res[2], int ifirst;

//    int zero;
//    int nedg;

//    int tsimp, tface, ssimp, sface;

//    int ncubes, first, last, k;
//    double u, v;
//    double p[2];
//    int type;
//    int half = 1;
//    int whole = 2;

//    int hn = 2;
//    int hm = 1;
//    int nsface_, nface_, nsoln_, nedges_;
//    int dims_ = 50;
//    int dime_ = 60;

//    double refval;
//    int i, j;

//    int ncvert_ = 4;
//    int nsimp_ = 2;

//    int numberOfCombinations = hc.combination(hn + 1, hm + 1);

//    //inicializing arrays
//    int storn_[hn + 1];
//    int storm_[hm + 1];
//    double cvert_[ncvert_][hn];
//    double vert[ncvert_][hn];
//    int bsvert_[hn + 1][hn];
//    int perm_[hn][nsimp_];
//    int comb_[numberOfCombinations][hm + 1];

//    //inicializing arrays dimensions
//    nsface_ = hc.mkcomb(&comb_[0][0], hn + 1, hm + 1);

//    int fnbr_[nsface_][nsface_];

//    int dimf_ = 5;

//    nsoln_ = -1;
//    double sol_[hn][dims_];
//    int solptr_[nsimp_][nsface_];
//    for (i = 0; i < nsimp_; i++) {
//        for (j = 0; j < nsface_; j++) solptr_[i][j] = 0;
//    } //TODO: Revisar como "solptr" eh modificada, os numero sao muito estranhos

//    int edges_[2][dime_];
//    for (i = 0; i < 2; i++) {
//        for (j = 0; j < dime_; j++) edges_[i][j] = -6;
//    } //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos

//    int smpedg_[nsimp_][2];
//    for (i = 0; i < nsimp_; i++) {
//        for (j = 0; j < 2; j++) smpedg_[i][j] = 0;
//    } //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos

//    //inicializing another arrays, it were globally defined in java
//    int facptr_[nsimp_][nsface_];
//    int face_[hm + 1][dimf_];

//    int dblev = 3;
//    printf("ifirst = %d\n", ifirst);
//    if (ifirst != 0) {
//        hc.mkcube(&cvert_[0][0], &bsvert_[0][0], &perm_[0][0], ncvert_, nsimp_, hn);
//        printf("Entering mkface():\n");
//        nface_ = hc.mkface(&face_[0][0], &facptr_[0][0], &fnbr_[0][0], dimf_, nsimp_, hn, hm, nsface_,
//                           &bsvert_[0][0], &comb_[0][0], &perm_[0][0], &storn_[0], &storm_[0]);
//        printf("After mkface(): nface = %d\n", nface_);

//    }

//    int exstfc[nface_];
//    for (i = 0; i < nface_; i++) exstfc[i] = 1;
//    int sptr_[nface_];

//    tsimp = 1;
//    tface = 3;

//    // End of initialization

//    // Curve2D proper
//    double foncub[hm][ncvert_];

//    // Work area
//    double g[hm][hm + 1];
//    double gx[hm];
//    double x[hm];
//    int wrki[hm];

//    // set the rectangle sizes and resolutions
//    double u0 = rect[0];
//    double u1 = rect[1];
//    double v0 = rect[2];
//    double v1 = rect[3];
//    int nu = res[0];
//    int nv = res[1];
//    double du = (u1 - u0) / nu;
//    double dv = (v1 - v0) / nv;

//    ncubes = nu * nv;

//    first = 1;
//    last = ncubes;

//    // initialize number of segments found
//    sn = 0;

//    //int i, j, k;

//    // loop over the rectangles
//    for (k = first; k <= last; k++) {
//        j = k % nv;
//        if (j == 0) {
//            j = nv;
//            i = k / nv;
//        } else {
//            i = (k / nv) + 1;
//        }
//        u = u0 + (i - 1) * du;
//        vert[0][0] = u + du;
//        vert[1][0] = u;
//        vert[2][0] = u + du;
//        vert[3][0] = u;

//        v = v0 + (j - 1) * dv;

//        vert[0][1] = v;
//        vert[1][1] = v;
//        vert[2][1] = v + dv;
//        vert[3][1] = v + dv;

//        // check whether all, half, or none of the square is inside

//        //lower right point
//        p[0] = u + du;
//        p[1] = v;

//        // a funcao inpdom foi criada a partir do arquivo bndry.F (localizada em phys/stone) do fortran!
//        //if (inpdom(&p[0]) == 0) goto lab200;
//        if (inpdom(&p[0]) != 0) {
//        // upper left point
//        p[0] = u;
//        p[1] = v + dv;
//        //if (inpdom(&p[0]) == 0) goto lab200;
//        if (inpdom(&p[0]) != 0) {

//        /* TODO: this works provided that the lower left corner is inside
//           when both the upper left and lower right corners are inside.
//           (e.g., rectangles oriented with the axes.)
//           upper right point */
//        p[0] = u + du;
//        p[1] = v + dv;

//        if (inpdom(&p[0]) == 1) {
//            type = whole;


//            ssimp = nsimp_;
//            sface = nface_;

//            printf("Inside if: sface = %d\n", sface);
//        } else {

//            type = half;
//            ssimp = tsimp;
//            sface = tface;
//            printf("Inside else: sface = %d\n", sface);
//        }

//        printf("Here: sface = %d\n", sface);

//        //set component 1
////        hc.putmf("VERT", &vert[0][0], ncvert_, n_);
//        zero = 0; // isso vai substituir o usso de uma variavel logica (false)

//        //TODO: lembrar descometar

//        //foncub[0][0] = f(vert[0][0], vert[0][1]);
//        RealVector u(2);
//        u(0) = vert[0][0];
//        u(1) = vert[0][1];

//        foncub[0][0] = hugoniot->HugoniotFunction(u);
//        //foncub[0][0] = f(vert[0][0], vert [0][1]);

//        refval = foncub[0][0];
//        for (int l = 1; l < 4; l++) {
//            if (l == 2 && type == half) goto lab90;
//              u(0) = vert[l][0];
//              u(1) = vert[l][1];

//              foncub[0][l] = hugoniot->HugoniotFunction(u);
//              //foncub[0][l] = f(vert[l][0], vert [l][1]);
//
//            if (refval * foncub[0][l] < 0.0) zero = 1;
//lab90:
//            ;
//        }

////        hc.putmf("FONCUB", &foncub[0][0], 4, 1);
//        //if (zero == 0) goto lab200;
//        if (zero != 0) {

//        //  solve for the current cube

//        // debug

//        //        if (dblev == 3) {
//        //            printf("CALLING CUBSOL WITH I, J=%d\n");
//        //        }
//        //DEBUG
//        //        status = hc.cubsol(&solptr_[0][0], &sol_[0][0], dims_, &sptr_[0], nsoln_,
//        //                &foncub[0][0], &vert[0][0], &exstfc[0], &face_[0][0], &facptr_[0][0],
//        //                hn, hm, ssimp, nsface_, sface, &storm_[0],
//        //                &storm_[0], &storm_[0]);
//        //

//        double u[hn][hm + 1]; //TODO Remove
//        double g[hm][hm+1];
//        double stormd[hm];
////        cout << "Valor de exstfc" << endl;
////        for (i = 0; i < nface_; i++) cout << exstfc[i] << " " << endl;

//        printf("Before cubsol(): sface = nface = %d\n", sface);
//        nsoln_ = hc.cubsol(&solptr_[0][0], &sol_[0][0], dims_, &sptr_[0], nsoln_, &foncub[0][0], &exstfc[0],
//                &face_[0][0], &facptr_[0][0], dimf_, &vert[0][0], ncvert_, hn, hm, ssimp, nsface_, sface, &u[0][0],
//                &g[0][0], &stormd[0], &storm_[0]);


//        //        hc.putmi("SOLPTR", &solptr_[0][0], nsimp_, nsface_);


//        // TODO: ver o numero certo de variavei de entrada e os tracinhos.

//        // debug
//        //        if ((dblev == 3) && (nsoln_ > 0)) {
////        hc.putmf("FONCUB", &foncub[0][0], ncvert_, m_);
////        hc.putmf("SOL", &sol_[0][0], n_, dims_);
////        hc.putmi("SPTR", &sptr_[0], nface_, 1);
////        hc.putmi("SOLPTR", &solptr_[0][0], nsimp_, nsface_);
////        //        }

//        // IMPROVE THE SOLUTION USING A ZERO-FINDER.

//        //        imprv(f, &sol_[0][0], dims_, sface, &sptr_[0], &face_[0][0], &vert[0][0]);
//        // (TODO: VER AS ENTRADAS CERTAS)

//        //MAKE THE LIST OF EDGE POINTERS

////        cout << "Valor de inteiros " << dime_ << " " << nedges_ << " " << nsimp_ << "  " << nsface_ << endl;



//        nedges_ = hc.mkedge(&edges_[0][0], dime_, nedges_, &smpedg_[0][0], &solptr_[0][0], &fnbr_[0][0], nsimp_, nsface_);
//        //for (int i = 0; i < 2; i++) for (int j = 0; j < dime_; j++) printf("edges_[%d][%d] = %d\n", i, j, edges_[i][j]);
////        printf("nedges_ = %d\n", nedges_);

////        printf("Estamos aqui\n");
////        cout << "Valor de nedges " << nedges_ << endl;

//        // debug
//        //        if ((dblev == 3) && (nsoln_ > 0)) {
////            hc.putmi("EDGES", &edges_[0][0], 2, nedges_);
////            hc.putmi("SMPEDG", &smpedg_[0][0], nsimp_, 2);
//        //        }

//        //STORE ALL PAIRS OF EDGES THAT WERE FOUND
//        //        cout << "Aqui sn= "<<sn<<endl;

//        if (nedges_ > 0) {
//            for (nedg = 0; nedg < nedges_; nedg++) {
////                printf("nedg = %d\n", nedg);
//                sn++;

////                segend[sn - 1][0][0] = sol_[0][edges_[0][nedg ]]; // X1
////                segend[sn - 1][0][1] = sol_[1][edges_[0][nedg ]]; // Y1
////                segend[sn - 1][1][0] = sol_[0][edges_[1][nedg ]]; // X2
////                segend[sn - 1][1][1] = sol_[1][edges_[1][nedg ]]; // Y2

////                cout << segend[sn - 1][0][0] << "    " << segend[sn - 1][0][1]<<";" << endl;

////                cout << segend[sn - 1][1][0] << "     " << segend[sn - 1][1][1] <<";"<< endl;

////                cout <<endl;

//                // Store the segments
////                printf("edges_[0][%d] = %d\n", nedg, edges_[0][nedg]);
////                printf("edges_[1][%d] = %d\n", nedg, edges_[1][nedg]);

//                RealVector p1(2), p2(2);
//                p1.component(0) = sol_[0][edges_[0][nedg ]]; //printf("1\n");
//                p1.component(1) = sol_[1][edges_[0][nedg ]]; //printf("2\n");

//                if (edges_[1][nedg ] < 0) return 0;
//                p2.component(0) = sol_[0][edges_[1][nedg ]]; //printf("3\n");
//                p2.component(1) = sol_[1][edges_[1][nedg ]]; //printf("4\n");

//                vrs.push_back(p1);
//                vrs.push_back(p2);

//                if (sn >= seglim) {
//                    return -1;
//                }
//            }

////            cout << "Valor de sn: " << sn << endl;
//        }
//
//
////lab200:;
//    }
// }
//}
//

//    }

//    printf("ENDS:   curv2d()\n\n");

    //contour2d();

    return 0;
}



int ContourMethodPure::contour2d(ImplicitFunction *impf, Boundary *boundary, double *rect, int *res, std::vector<RealVector> &vrs) {

//TODO: int sn must be eliminated (because it is == vrs.size()).
// TODO: Get rid of seglim & fdummy, ifirst (use ctor instead), and rect will become Domain*.
//    printf("BEGINS: vect2d()\n");
//    HyperCube hc;
//    deallocate_arrays();
    allocate_arrays();

    vrs.clear();

    int zero;
    int nedg;

//    int tsimp, tface, ssimp, sface;

    int ssimp, sface;

    double u, v;
    double p[2];
    int type;
    int half = 1;
    int whole = 2;

//    int hn = 2;
//    int hm = 1;
//    int nsface_, nface_, nsoln_, nedges_;
//    int dims_ = 50;
//    int dime_ = 60;
//    int dimf_ = 5;
//    int ncvert_ = 4;
//    int nsimp_ = 2;

//    int numberOfCombinations = hc.combination(hn + 1, hm + 1); printf("numberOfCombinations = %d\n", numberOfCombinations);

    double refval;
    int i, j;

    //inicializing arrays
//    int storn_[hn + 1];
//    int storm_[hm + 1];
//    double cvert_[ncvert_][hn];
//    double vert[ncvert_][hn];
//    int bsvert_[hn + 1][hn];
//    int perm_[hn][nsimp_];
//    int comb_[numberOfCombinations][hm + 1];

    //inicializing arrays dimensions
//    nsface_ = hc.mkcomb(comb_, hn + 1, hm + 1); printf("nsface_ = %d\n", nsface_);

//    int fnbr_[nsface_][nsface_];

    //nsoln_ = -1;
//    double sol_[hn][dims_];
//    int solptr_[nsimp_][nsface_];
//    for (i = 0; i < nsimp_; i++) {
//        //for (j = 0; j < nsface_; j++) solptr_[i][j] = -150;
//    } //TODO: Revisar como "solptr" eh modificada, os numero sao muito estranhos

//    int edges_[2][dime_];
//    for (i = 0; i < 2; i++) {
//        //for (j = 0; j < dime_; j++) edges_[i][j] = -6;
//    } //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos

//    int smpedg_[nsimp_][2];
//    for (i = 0; i < nsimp_; i++) {
//        //for (j = 0; j < 2; j++) smpedg_[i][j] = 0;
//    } //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos

//    //inicializing another arrays, it were globally defined in java
//    int facptr_[nsimp_][nsface_];
//    int face_[hm + 1][dimf_];

//    int dblev = 3;
//    printf("ifirst = %d\n", ifirst);
//    if (ifirst != 0) {
//        hc.mkcube(cvert_, bsvert_, perm_, ncvert_, nsimp_, hn);
//        printf("Entering mkface():\n");
//        nface_ = hc.mkface(&face_[0][0], &facptr_[0][0], fnbr_, dimf_, nsimp_, hn, hm, nsface_,
//                           bsvert_, comb_, perm_, storn_, storm_);
//        printf("After mkface(): nface = %d\n", nface_);

//    }

//    int exstfc[nface_];
//    //for (i = 0; i < nface_; i++) exstfc[i] = 1;
//    int sptr_[nface_];

//    tsimp = 1;
//    tface = 3;

    // End of initialization

    // Curve2D proper
    double foncub[hm][ncvert_];

//    // Work area
//    double g[hm][hm + 1];
//    double gx[hm];
//    double x[hm];
//    int wrki[hm];

    // set the rectangle sizes and resolutions
    double u0 = rect[0];
    double u1 = rect[1];
    double v0 = rect[2];
    double v1 = rect[3];
    int    nu = res[0];
    int    nv = res[1];
    double du = (u1 - u0) / nu;
    double dv = (v1 - v0) / nv;

//    printf("nu = %d, nv = %d\n", nu, nv);
//    printf("du = %f, dv = %f\n", du, dv);

    for (i = 0; i < nface_; i++) exstfc[i] = 1;

    // initialize number of segments found
//    sn = 0;

    // loop over the rectangles
    /* Recall that the order of vertex index on the rectangles is given by hcube and it is not trivial,
       the order is:
                       3-vertex +---+ 2-vertex
                                |\  |
                                | \ |
                                |  \|
                       1-vertex +---+ 0-vertex
    */
    for (i = 0; i < nu; i++) {
        u = u0 + i * du;
        vert[0*hn + 0] = u + du;
        vert[1*hn + 0] = u;
        vert[2*hn + 0] = u + du;
        vert[3*hn + 0] = u;


//        cout << "X: "<<  vert[0*hn + 0]<< " "<<vert[1*hn + 0]<< "  "<<vert[2*hn + 0]  <<"  "<< vert[3*hn + 0] <<endl;

        for (j = 0; j < nv; j++) {
            v = v0 + j * dv;

            vert[0*hn + 1] = v;
            vert[1*hn + 1] = v;
            vert[2*hn + 1] = v + dv;
            vert[3*hn + 1] = v + dv;

//            cout << "Y: "<< vert[0 * hn + 1] << " " << vert[1 * hn + 1] << "  " << vert[2 * hn + 1] << "  " << vert[3 * hn + 1] << endl;



            // check whether all, half, or none of the square is inside

            //lower right point
            p[0] = u + du;
            p[1] = v;

            // a funcao inpdom foi criada a partir do arquivo bndry.F (localizada em phys/stone) do fortran!
            //if (inpdom(&p[0]) == 0) goto lab200;
//            if (inpdom(&p[0]) != 0) {
            if (boundary->inside(&p[0])) {
                // upper left point
                p[0] = u;
                p[1] = v + dv;
                //if (inpdom(&p[0]) == 0) goto lab200;
//                if (inpdom(&p[0]) != 0) {
                if (boundary->inside(&p[0])) {

                    /* TODO: this works provided that the lower left corner is inside
                       when both the upper left and lower right corners are inside.
                       (e.g., rectangles oriented with the axes.)
                       upper right point */
                    p[0] = u + du;
                    p[1] = v + dv;

//                    if (inpdom(&p[0]) == 1) {
                    if (boundary->inside(&p[0])) {
                        type = whole;
                        ssimp = nsimp_;
                        sface = nface_;
//                        printf("Inside if: sface = %d\n", sface);
                    } else {
                        type = half;
                        ssimp = tsimp;
                        sface = tface;
//                        printf("Inside else: sface = %d\n", sface);
                    }

//                    printf("Here: sface = %d\n", sface);

                    //set component 1
                    zero = 0; // isso vai substituir o usso de uma variavel logica (false) TODO : --> colocar como booleano (zero --> root)

                    // foncub is filled here
//                    int status = impf->function_on_square(&foncub[0][0], i, j, type);
                    int status = impf->function_on_square(&foncub[0][0], i, j);
//                    printf("FONCUB = %f, %f, %f, %f\n",foncub[0][0],foncub[0][1],foncub[0][2],foncub[0][3]);


                    if (status != 0) {
                        refval = foncub[0][0];
                        for (int l = 1; l < 4; l++) {
                            if (l == 2 && type == half) {
                                //    goto lab90;
                                } else {
                                    if (refval * foncub[0][l] < 0.0) zero = 1;
                                //    lab90:
                                //    ;
                                }
                            }
                        }

//                    printf("Here: zero = %d\n", zero);

                    if (zero != 0) {
                        //  solve for the current cube

                        // debug
                        //        if (dblev == 3) {
                        //            printf("CALLING CUBSOL WITH I, J=%d\n");
                        //        }
                        // DEBUG

                        double u[hn][hm + 1]; //TODO Remove
                        double g[hm][hm+1];
                        double stormd[hm];

//                        printf("Before cubsol(): sface = nface = %d\n", sface);
                        nsoln_ = hc.cubsol(solptr_, sol_, dims_, sptr_, nsoln_, &foncub[0][0],
                                           exstfc, face_, facptr_, dimf_, vert,
                                           ncvert_, hn, hm, ssimp, nsface_, sface, &u[0][0], &g[0][0],
                                           stormd, storm_);

                        // TODO: ver o numero certo de variavei de entrada e os tracinhos.

                        // debug
                        //        if ((dblev == 3) && (nsoln_ > 0)) {
                        //            hc.putmf("FONCUB", &foncub[0][0], ncvert_, m_);
                        //            hc.putmf("SOL", &sol_[0][0], n_, dims_);
                        //            hc.putmi("SPTR", &sptr_[0], nface_, 1);
                        //            hc.putmi("SOLPTR", &solptr_[0][0], nsimp_, nsface_);
                        //        }
                        // DEBUG

                        // IMPROVE THE SOLUTION USING A ZERO-FINDER.
                        //        imprv(f, &sol_[0][0], dims_, sface, &sptr_[0], &face_[0][0], &vert[0][0]);
                        // (TODO: VER AS ENTRADAS CERTAS)

                        // DESTROY BELOW //
                        if (impf->improvable()) {
                            for (int ii = 0; ii < sface; ii++){

                                if (sptr_[ii] == -1) continue;
                                int sp = sptr_[ii];

                                RealVector p0newton(2), p1newton(2), p_init_newton(2), p_improved_newton(2);

                                p0newton.component(0) = vert[face_[0*dimf_ + ii]*hn + 0]; // p1(1) = vert(1,face(1,i)+1)
                                p0newton.component(1) = vert[face_[0*dimf_ + ii]*hn + 1]; // p1(2) = vert(2,face(1,i)+1)

                                p1newton.component(0) = vert[face_[1*dimf_ + ii]*hn + 0]; // p2(1) = vert(1,face(2,i)+1)
                                p1newton.component(1) = vert[face_[1*dimf_ + ii]*hn + 1]; // p2(2) = vert(2,face(2,i)+1)

                                // To initialize:
                                //
                                for (int jj = 0; jj < 2; jj++) p_init_newton.component(jj) = sol_[jj*dims_ + sp];

                                Newton_Improvement *newton_improver = new Newton_Improvement(impf);
                                int newton_info = newton_improver->newton(p0newton, p1newton, p_init_newton, p_improved_newton);
//                                if (newton_info == NEWTON_IMPROVEMENT_ERROR){
//                                    printf("        dimf_ = %d, face_[%d] = %d\n", dimf_, 0*dimf_ + ii, face_[0*dimf_ + ii]);
//                                }

//printf("p0 = (%lf,%lf), p1 = (%lf,%lf), p_init = (%lf,%lf)\n",
// p0newton.component(0), p0newton.component(1), p1newton.component(0),
// p1newton.component(1), p_init_newton.component(0), p_init_newton.component(1));

                                sol_[0*dims_ + sp] = p_improved_newton.component(0); // sol(1,sp) = v(1)
                                sol_[1*dims_ + sp] = p_improved_newton.component(1); // sol(2,sp) = v(2)

                            }
                        }
                        // DESTROY ABOVE //


                        //MAKE THE LIST OF EDGE POINTERS
                        int N_EDGES = hc.mkedge(edges_, dime_, nedges_, smpedg_, solptr_,
                                            fnbr_, ssimp, nsface_); // TODO: Em Fortran esta nsimp_

                        nedges_ = N_EDGES;
//                        printf("nedges_ = %d\n", nedges_);

                        // debug
                        //        if ((dblev == 3) && (nsoln_ > 0)) {
                        //            hc.putmi("EDGES", &edges_[0][0], 2, nedges_);
                        //            hc.putmi("SMPEDG", &smpedg_[0][0], nsimp_, 2);
                        //        }
                        // DEBUG

//                               printf("nedges = %d\n", nedges_);

                        //STORE ALL PAIRS OF EDGES THAT WERE FOUND
                        if (nedges_ > 0) {
                           for (nedg = 0; nedg < nedges_; nedg++) {
//                               printf("nedg = %d\n", nedg);
//                               sn++;

                               // Store the segments
//                               printf("edges_[0][%d] = %d\n", nedg, edges_[0][nedg]);
//                               printf("edges_[1][%d] = %d\n", nedg, edges_[1][nedg]);

                               RealVector p1(2), p2(2);
                               p1.component(0) = sol_[0*dims_ + edges_[0*dime_ + nedg ]]; //printf("1\n");
                               p1.component(1) = sol_[1*dims_ + edges_[0*dime_ + nedg ]]; //printf("2\n");

// TODO: Pablo comentou a seguinte linha em 11 Janeiro 2012
//                               if (edges_[1*dime_ + nedg ] < 0) return 0;
                               p2.component(0) = sol_[0*dims_ + edges_[1*dime_ + nedg ]]; //printf("3\n");
                               p2.component(1) = sol_[1*dims_ + edges_[1*dime_ + nedg ]]; //printf("4\n");

                               vrs.push_back(p1);
                               vrs.push_back(p2);

//                               if (sn >= seglim) {
//                                   return -1;
//                               }
                            }
                        }
                    }
                }
            }
        }
    }
//    printf("ENDS:   vect2d()\n\n");
    return 0;
}


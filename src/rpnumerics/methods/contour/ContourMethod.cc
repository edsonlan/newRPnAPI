/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) HugoniotContourMethod.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "ContourMethod.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

// Static variables are defined here:
//
bool ContourMethod::is_first = true;

HyperCube ContourMethod::hc;

int ContourMethod::hn = 2;
int ContourMethod::hm = 1;

int ContourMethod::nsface_ = 3;
int ContourMethod::nface_ = 5;
int ContourMethod::nsoln_ = -1;
int ContourMethod::nedges_;

int ContourMethod::dims_ = 6; // Estava 50, Fortran diz 6 (inicialmente 4)
int ContourMethod::dime_ = 6; // Estava 60, Fortran diz 6 (inicialmente 2)
int ContourMethod::dimf_ = 5;
int ContourMethod::ncvert_ = 4;
int ContourMethod::nsimp_ = 2;

int ContourMethod::numberOfCombinations = 3;

int * ContourMethod::storn_;
int * ContourMethod::storm_;

double * ContourMethod::cvert_;
double * ContourMethod::vert;
int * ContourMethod::bsvert_;
int * ContourMethod::perm_;
int * ContourMethod::comb_;
int * ContourMethod::fnbr_;
int * ContourMethod::facptr_;
int * ContourMethod::face_;
double * ContourMethod::sol_;
int * ContourMethod::solptr_;
int * ContourMethod::edges_;
int * ContourMethod::smpedg_;
int * ContourMethod::exstfc;
int * ContourMethod::sptr_;
int * ContourMethod::gamb;

int ContourMethod::tsimp = 1;
int ContourMethod::tface = 3;

Matrix<int>  ContourMethod::number_chains;
Matrix<std::vector <std::vector <int> > > ContourMethod::chain_edges;
Matrix<std::vector <std::vector <RealVector> > > ContourMethod::chains;
// Matrix<bool> ContourMethod::iplus, ContourMethod::iminus; 
// Matrix<bool> ContourMethod::jplus, ContourMethod::jminus;

std::vector <std::vector <int> > ContourMethod::chain_list;


void ContourMethod::allocate_arrays(void){
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

        exstfc = new int[nface_]; for (int i = 0; i < nface_; i++) exstfc[i] = 1;
        sptr_ = new int[nface_];

        // This is a useful GAMBIARRA
        gamb = new int[5];


        is_first = false;

        printf("++++++++++++++++ ContourMethod: REMEMBER TO INVOKE deallocate_arrays() AT QUIT-TIME!!!\n++++++++++++++++ DON\'T SAY I DIDN\'T WARN YOU!!!\n");
    }

    return;
}


void ContourMethod::deallocate_arrays(void){
    if (!is_first){
        delete [] gamb;
        delete [] sptr_;
        delete [] exstfc;

        delete [] smpedg_;
        delete [] edges_;
        delete [] solptr_;
        delete [] sol_;

        delete [] face_;
        delete [] facptr_;

        delete [] fnbr_;

        delete [] comb_;
        delete [] perm_;
        delete [] bsvert_;
        delete [] vert;
        delete [] cvert_;

        delete [] storm_;
        delete [] storn_;

        is_first = true;
    }

    return;
}


int ContourMethod::contour2d(ImplicitFunction *impf, std::vector<RealVector> &vrs) {
    std::vector< std::deque <RealVector> > curves;
    std::vector < bool > is_circular;
    int method = SEGMENTATION_METHOD;

    return contour2d(impf, vrs, curves, is_circular, method);
}

int ContourMethod::contour2d(ImplicitFunction *impf, std::vector< std::deque <RealVector> > &curves, std::vector <bool> &is_circular) {
    std::vector<RealVector> vrs;
    int method = CONTINUATION_METHOD;

    return contour2d(impf, vrs, curves, is_circular, method);
}

int ContourMethod::contour2d(ImplicitFunction *impf, int min_row, int max_row, int min_col, int max_col,
                             std::vector<RealVector> &vrs,
                             std::vector< std::deque <RealVector> > &curves, std::vector <bool> &is_circular,
                             const int method) {

    GridValues *gv = impf->grid_value();

//    printf("BEGINS: Contour2D()\n");

//    deallocate_arrays();
    allocate_arrays();

//    vrs.clear();
//    curves.clear();
//    is_circular.clear();

    int zero;
    int nedg;

    int ssimp, sface;

    double refval;

    // The curve is tested in small squares with values on foncub
    double foncub[hm][ncvert_];

    int rows = gv->grid.rows() - 1;
    int cols = gv->grid.cols() - 1;


// There are two methods, when the search algorithm is concatenated, some new structures are need.
// Here the structures initialize. [number_chains(i,j) = 0] is used as a flag below.
if ( method == CONTINUATION_METHOD ) {
    number_chains.resize(rows,cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            number_chains(i,j) = 0;
        }
    }

    chains.resize(rows,cols);
    chain_edges.resize(rows,cols);
//    iplus.resize(rows,cols);
//    iminus.resize(rows,cols);
//    jplus.resize(rows,cols);
//    jminus.resize(rows,cols);
    chain_list.clear();
}

    std::vector<int> coordinates(2);

    /* Recall that the order of vertex index on the rectangles is given by hcube;
       it is not trivial, the order is:
                       3-vertex +---+ 2-vertex
                                |\  |
                                | \ |
                                |  \|
                       1-vertex +---+ 0-vertex

       The order of the faces is also given by hcube:
                                4-face
                                +----+
                                |\   |
                        2-face  | \  |  3-face
                                |  \ |
                                +----+
                                0-face
       The diagonal connectic 0-vertex to 3-vertex is 1-face.
    */

//    std::ofstream off("cell.txt", std::ifstream::out);

    // loop over the rectangles
    for (int i = min_row; i < max_row; i++) { //std::cout << "i = " << i << std::endl;
        for (int j = min_col; j < max_col; j++) { //std::cout << "    j = " << j << std::endl;
            vert[0*hn + 0] = gv->grid(i + 1, j + 0).component(0); // u + du;
            vert[1*hn + 0] = gv->grid(i + 0, j + 0).component(0); // u;
            vert[2*hn + 0] = gv->grid(i + 1, j + 1).component(0); // u + du;
            vert[3*hn + 0] = gv->grid(i + 0, j + 1).component(0); // u;

            vert[0*hn + 1] = gv->grid(i + 1, j + 0).component(1); // v;
            vert[1*hn + 1] = gv->grid(i + 0, j + 0).component(1); // v;
            vert[2*hn + 1] = gv->grid(i + 1, j + 1).component(1); // v + dv;
            vert[3*hn + 1] = gv->grid(i + 0, j + 1).component(1); // v + dv;

            //std::cout << "    After verts" << std::endl;

//            off << "    Cell " << i << ", " << j << " is: " << gv->cell_type(i, j) << std::endl;

            if (gv->cell_type(i, j) != CELL_IS_INVALID) {
                    if (gv->cell_type(i, j) == CELL_IS_SQUARE){
                        ssimp = nsimp_;
                        sface = nface_;
                    } 
                    else if (gv->cell_type(i, j) == CELL_IS_TRIANGLE) {
                        ssimp = tsimp;
                        sface = tface;
                    }
                    else continue;

                    // The posibility of a zero on a cell is given by this...
                    zero = 0; // isso vai substituir o usso de uma variavel logica (false) TODO : --> colocar como booleano (zero --> root)

                    // foncub is filled here
                    //std::cout << "    before status..." << std::endl;
                    int status = impf->function_on_square(&foncub[0][0], i, j);
                    //std::cout << "    after status = " << status << std::endl;
                    
                    if (status != 0) {
                        refval = foncub[0][0];
                        for (int l = 1; l < 4; l++) {
                            if (l == 2 && gv->cell_type(i, j) == CELL_IS_TRIANGLE) {
                            }
                            else {
                                if (refval * foncub[0][l] < 0.0) zero = 1;
                            }
                        }
                    }

                    //std::cout << "    zero = " << zero << std::endl;

                    if (zero != 0) {
                        //  solve for the current cube
                        double u[hn][hm + 1]; //TODO Remove
                        double g[hm][hm+1];
                        double stormd[hm];

                        nsoln_ = hc.cubsol(solptr_, sol_, dims_, sptr_, nsoln_, &foncub[0][0],
                                           exstfc, face_, facptr_, dimf_, vert,
                                           ncvert_, hn, hm, ssimp, nsface_, sface, &u[0][0], &g[0][0],
                                           stormd, storm_);

                        // IMPROVE THE SOLUTION USING A ZERO-FINDER.
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

                                sol_[0*dims_ + sp] = p_improved_newton.component(0); // sol(1,sp) = v(1)
                                sol_[1*dims_ + sp] = p_improved_newton.component(1); // sol(2,sp) = v(2)

                            }
                        }

                        // MAKE THE LIST OF EDGE POINTERS
                        int N_EDGES = hc.mkedge(edges_, dime_, nedges_, smpedg_, solptr_,
                                            fnbr_, ssimp, nsface_); // TODO: Em Fortran esta nsimp_

                        nedges_ = N_EDGES;

                        //STORE ALL PAIRS OF EDGES THAT WERE FOUND
                        if (nedges_ > 0) {
                            // A Topological Sort is called for ordering the segments in each cell.
                            //
                            if ( method == CONTINUATION_METHOD ) {
                                if( !topological_sort(i,j) ) continue;
                                coordinates[0] = i; coordinates[1] = j;
                                // Store the coodinates of the GridValues that have solution, so after
                                // the sort, a list of order n is used and not the whole grid n^2.
                                chain_list.push_back(coordinates);
//                                cout << "Coordinates: " << coordinates[0] << ", " << coordinates[1] << endl;
//                                cout << "Number of chains: " << number_chains(i,j) << endl;
//                                cout << "The first point:  " << chains(i,j)[0][0] << endl;
                            }

                            // The segments are store in pairs
                            //
                            if ( method == SEGMENTATION_METHOD ) {
                                for (nedg = 0; nedg < nedges_; nedg++) {

                                   // Store the segments
                                   RealVector p1(2), p2(2);
                                   p1.component(0) = sol_[0*dims_ + edges_[0*dime_ + nedg ]];
                                   p1.component(1) = sol_[1*dims_ + edges_[0*dime_ + nedg ]];

                                   // TODO: Pablo comentou a seguinte linha em 11 Janeiro 2012
                                   //       if (edges_[1*dime_ + nedg ] < 0) return 0;
                                   p2.component(0) = sol_[0*dims_ + edges_[1*dime_ + nedg ]];
                                   p2.component(1) = sol_[1*dims_ + edges_[1*dime_ + nedg ]];

                                   vrs.push_back(p1);
                                   vrs.push_back(p2);
                                }
                            }
                        }
                    }

                    // 
            }
        }
    }

//    off.close();

    if ( method == SEGMENTATION_METHOD ) {
//        printf("ENDS:     Contour2D()/For Segments\n\n");
        return 0;
    }

// TEST
//   cout << "Cells: ";
//   for (int i = 0; i < chain_list.size(); i++) {
//       cout << number_chains(chain_list[i][0],chain_list[i][1]) << " ";
//   }
//   cout << endl;
// End TEST

    // The following integers keep the path of the concatenation
    int prev, middle, next;

    curves.clear();
    is_circular.clear();

    // The loop is given over the active cell grid stored in [chain_list]
    for (int i = 0; i < chain_list.size(); /*For two segments, the cell is used twice*/) {
        int i_index = chain_list[i][0]; int i_start = i_index;
        int j_index = chain_list[i][1]; int j_start = j_index;

        // Here is the increment of [int i], when needed
        if ( number_chains(i_index, j_index) < 2 ) i++;
        if ( number_chains(i_index, j_index) == 0 ) continue;

        // Reused for each branch, is a <deque> because it can be filled by both sides.
        std::deque <RealVector> branch;
        bool circular_is = false;
        int chain_size = chains(i_index,j_index)[0].size();

        // [prev] store the direction not followed, [middle] the direction to follow; [next] will
        // store the entering direction in the next cell.
        prev   = chain_edges(i_index,j_index)[0][0];
        middle = chain_edges(i_index,j_index)[0][1];

        // Store the points in the [branch] and eliminate the points of the list 
        for (int j = 0; j < chain_size; j++) branch.push_back( chains(i_index,j_index)[0][j] );
        chains(i_index,j_index)[0].erase( chains(i_index,j_index)[0].begin());
        chain_edges(i_index,j_index)[0].erase( chain_edges(i_index,j_index)[0].begin());

        do {
            restore:

            // [middle] also works as a flag. If it is used once then it is setted as -1, there are
            //          no previous direction, then returning to the same cell end the procedure.
            if (middle == -1) break;
            number_chains(i_index, j_index)--;

            // [middle] gives the direction to follow and the index of the next point
            switch ( middle ) {
                case 0: {
                    j_index--;
                    next = 4;
                }
                break;

                case 2: {
                    i_index--;
                    next = 3;
                }
                break;

                case 3: {
                    i_index++;
                    next = 2;
                }
                break;

                case 4: {
                    j_index++;
                    next = 0;
                }
                break;

                default:
                break;
            }
            // Preventing to leave the domain.
            if ( (i_index < 0) || (j_index < 0) || (i_index >= rows) || (j_index >= cols) || middle == 1 ) {
                middle = prev;
                i_index = i_start;
                j_index = j_start;
                prev = -1;
//                /* DEBUG... */
//                cout << " :: FIZ __RESTORE__" << endl;
//                cout << "index(" << next << "): " << i_index << ", " << j_index << endl;
//                cout << "i: " << i_index << ", " << i_start << endl;
//                cout << "j: " << j_index << ", " << j_start << endl;
//                /* ... END DEBUG */
                goto restore;
            }

            // Some restriction may cut the continuation of the curve, the following IF saves this issues.
            if( chain_edges(i_index,j_index).size() == 0 ) break;
            // Store the points in the [branch] and eliminate the points of the list.
            // The chain in (i, j) must be identified with [next] the direction of push_back or
            // push_front is denoted by the use of [prev = -1], else it will be the middle direction.
            for (int k = 0; k < chains(i_index,j_index)[0].size(); k++) {
                if ( chain_edges(i_index,j_index)[k][0] == next ) {
                    chain_size = chains(i_index,j_index)[k].size();
                    if (prev == -1)
                        for (int j = 1; j < chain_size; j++) branch.push_front( chains(i_index,j_index)[k][j] );
                    else
                        for (int j = 1; j < chain_size; j++) branch.push_back( chains(i_index,j_index)[k][j] );
                    middle = chain_edges(i_index,j_index)[k][1];
                    chains(i_index,j_index).erase( chains(i_index,j_index).begin() + k );
                    chain_edges(i_index,j_index).erase( chain_edges(i_index,j_index).begin() + k);
                    if ( (i_index == i_start) && (j_index == j_start) && (next == prev) ) {
                        // TODO: The repeated point it is not deleted! It is convenient to keep it.
                        // branch.erase( branch.begin() );
                        circular_is = true;
                    }
                    break;
                }
                if ( chain_edges(i_index,j_index)[k][1] == next ) {
                    chain_size = chains(i_index,j_index)[k].size();
                    if (prev == -1) 
                        for (int j = chain_size-2; j >= 0; j--) branch.push_front( chains(i_index,j_index)[k][j] );
                    else            
                        for (int j = chain_size-2; j >= 0; j--) branch.push_back( chains(i_index,j_index)[k][j] );
                    middle = chain_edges(i_index,j_index)[k][0];
                    chains(i_index,j_index).erase( chains(i_index,j_index).begin() + k );
                    chain_edges(i_index,j_index).erase( chain_edges(i_index,j_index).begin() + k);
                    if ( (i_index == i_start) && (j_index == j_start) && (next == prev) ) {
                        // TODO: The repeated point it is not deleted! It is convenient to keep it.
                        // branch.erase( branch.begin() );
                        circular_is = true;
                    }
                    break;
                }
            }

        } while ( number_chains(i_index, j_index) > 0 );

        curves.push_back(branch);
        is_circular.push_back(circular_is);
    }

// For DEBUGING, the following lines can be included instead of the two ereasing statements.
///* DEBUG...
//cout << "Antes de deletar tenho: "<< endl;
//for(int ii = 0; ii < chains(i_index,j_index).size(); ii++){
//    for (int jj = 0; jj < chains(i_index,j_index)[ii].size(); jj++){
//        cout << " " << chains(i_index,j_index)[ii][jj];
//    }
//    cout << endl;
//}
//... END DEBUG */
//                    chains(i_index,j_index).erase( chains(i_index,j_index).begin() + k );
//                    chain_edges(i_index,j_index).erase( chain_edges(i_index,j_index).begin() + k);
///* DEBUG...
//cout << "Depois de deletar tenho: "<< endl;
//for(int ii = 0; ii < chains(i_index,j_index).size(); ii++){
//    for (int jj = 0; jj < chains(i_index,j_index)[ii].size(); jj++){
//        cout << " " << chains(i_index,j_index)[ii][jj];
//    }
//    cout << endl;
//}
//... END DEBUG */

///*DEBUG: Impressions*/
//    cout << "Finalizei com " << curves.size() << " curvas :: ::   ";
//    for (int k = 0; k < curves.size(); k++) cout << "Tamanho de " << k << " es: " << curves[k].size() << " y es circular: " << is_circular[k] << "; ";
//    cout << endl;

//    for (int k = 0; k < curves.size(); k++) {
//        cout << "Curva " << k << ":";
//        for (int l = 0; l < curves[k].size(); l++) cout << " " << curves[k][l];
//        cout << endl;
//    }
///*END DEBUG...*/

    printf("ENDS:     Contour2D()/Continuous\n\n");

    return 0;
}


/**
 *  This method is the key to the segments ordering finded by the search algorithm.
 *  There are some tricks :
 *     [number_chains] stores the two possible cases, 1 segment or 1 chain or two chains given
 *                     by a segment each.
 *     [chains] has the dimension given by the number of chains, each chain is given separetally.
 *     [chain_edges] uses the useful GAMBIARRA which for each initial and final point gives the
 *                   index of the cell face reached.
**/

int ContourMethod::topological_sort(int i, int j) {
        for (int ii = 0; ii < 5; ii++) {
            gamb[ii] = 9;
            for (int jj = 0; jj < 5; jj++) {
                if (ii == sptr_[jj]) { 
                    gamb[ii] = jj;
                    break;
                }
            }
        }
    /* Aqui entra o topological sort */
    if (nedges_ == 2) {
        if (edges_[dime_ + 1] == 3) {
            // TODO: Sabemos exatamente quem eh quem
            number_chains(i,j) = 2;

            chains(i,j).resize(2);
            chains(i,j)[0].resize(2);
            chains(i,j)[1].resize(2);
            chains(i,j)[0][0].resize(2);
            chains(i,j)[0][1].resize(2);
            chains(i,j)[1][0].resize(2);
            chains(i,j)[1][1].resize(2);

            for (int k = 0; k < 2; k++) chains(i,j)[0][0][k] = sol_[k*dims_ + 0];
            for (int k = 0; k < 2; k++) chains(i,j)[0][1][k] = sol_[k*dims_ + 1];

            for (int k = 0; k < 2; k++) chains(i,j)[1][0][k] = sol_[k*dims_ + 2];
            for (int k = 0; k < 2; k++) chains(i,j)[1][1][k] = sol_[k*dims_ + 3];

            chain_edges(i,j).resize(2);
            chain_edges(i,j)[0].resize(2);
            chain_edges(i,j)[1].resize(2);
            // iminus e jminus
            chain_edges(i,j)[0][0] = gamb[0]; chain_edges(i,j)[0][1] = gamb[1];
            // jplus e iplus
            chain_edges(i,j)[1][0] = gamb[2]; chain_edges(i,j)[1][1] = gamb[3];

            // iminus(i,j) = iplus(i,j) = jminus(i,j) = jplus(i,j) = true;
        }
        else if ( (edges_[1] + edges_[dime_ + 1]) == 3) {
            //cout << "Ordem: 0 1 2: 3 = " << edges_[1] << " + " << edges_[dime_ + 1] << endl;
            // TODO: Podemos inferir quem eh quem
            number_chains(i,j) = 1;
            chains(i,j).resize(1);
            chains(i,j)[0].resize(3);

            chains(i,j)[0][0].resize(2);
            chains(i,j)[0][1].resize(2);
            chains(i,j)[0][2].resize(2);

            for (int k = 0; k < 2; k++) chains(i,j)[0][0][k] = sol_[k*dims_ + 0];
            for (int k = 0; k < 2; k++) chains(i,j)[0][1][k] = sol_[k*dims_ + 1];
            for (int k = 0; k < 2; k++) chains(i,j)[0][2][k] = sol_[k*dims_ + 2];

            chain_edges(i,j).resize(1);
            chain_edges(i,j)[0].resize(2);
            // chain = {0, 1, 2}
            chain_edges(i,j)[0][0] = gamb[0]; chain_edges(i,j)[0][1] = gamb[2];
        }
        else {
            //cout << "Ordem: 1 0 2: 2 = " << edges_[1] << " + " << edges_[dime_ + 1] << endl;
            // TODO: Sabemos exatamente quem eh quem
            number_chains(i,j) = 1;
            chains(i,j).resize(1);
            chains(i,j)[0].resize(3);

            chains(i,j)[0][0].resize(2);
            chains(i,j)[0][1].resize(2);
            chains(i,j)[0][2].resize(2);

            for (int k = 0; k < 2; k++) chains(i,j)[0][0][k] = sol_[k*dims_ + 1];
            for (int k = 0; k < 2; k++) chains(i,j)[0][1][k] = sol_[k*dims_ + 0];
            for (int k = 0; k < 2; k++) chains(i,j)[0][2][k] = sol_[k*dims_ + 2];

            chain_edges(i,j).resize(1);
            chain_edges(i,j)[0].resize(2);
            // chain = {1, 0, 2}
            chain_edges(i,j)[0][0] = gamb[1]; chain_edges(i,j)[0][1] = gamb[2];
        }
    }
    else {
        // TODO: Podemos inferir quem eh quem
        number_chains(i,j) = 1;
        chains(i,j).resize(1);
        chains(i,j)[0].resize(2);

        chains(i,j)[0][0].resize(2);
        chains(i,j)[0][1].resize(2);

        for (int k = 0; k < 2; k++) chains(i,j)[0][0][k] = sol_[k*dims_ + 0];
        for (int k = 0; k < 2; k++) chains(i,j)[0][1][k] = sol_[k*dims_ + 1];

        chain_edges(i,j).resize(1);
        chain_edges(i,j)[0].resize(2);
        // chain = {0, 1}
        chain_edges(i,j)[0][0] = gamb[0]; chain_edges(i,j)[0][1] = gamb[1];
    }
    /* Aqui termina o tsort */

/* DEBUG...
    cout << "Impressoes para (" << i << ", " << j << "): " << number_chains(i,j) << " Chains" << endl;
    cout << "sptr_:";
    for (int j = 0; j < 5; j++) cout << " " << sptr_[j];
    cout << endl << "gamb: ";
    for (int j = 0; j < 5; j++) cout << " " << gamb[j];
    cout << endl;

    for (int ii = 0; ii < number_chains(i,j); ii++) {
        cout << " cadeia(" << ii+1 << "): ";
        for (int jj = 0; jj < chains(i,j)[ii].size(); jj++) cout << chains(i,j)[ii][jj] << " ";
        cout << endl;
        cout << " finales: " << chain_edges(i,j)[ii][0] << " -- " << chain_edges(i,j)[ii][1] << endl;
    }
    cout << endl;
... DEBUG*/

    return 1;
}

int ContourMethod::contour2d(ImplicitFunction *impf, 
                         std::vector<RealVector> &vrs,
                         std::vector< std::deque <RealVector> > &curves,
                         std::vector <bool> &is_circular,
                         const int method){
    GridValues *gv = impf->grid_value();

//    printf("BEGINS: Contour2D()\n");

//    deallocate_arrays();
    allocate_arrays();

    vrs.clear();
    curves.clear();
    is_circular.clear();

    int zero;
    int nedg;

    int ssimp, sface;

    double refval;

    // The curve is tested in small squares with values on foncub
    double foncub[hm][ncvert_];

    int rows = gv->grid.rows() - 1;
    int cols = gv->grid.cols() - 1;


// There are two methods, when the search algorithm is concatenated, some new structures are need.
// Here the structures initialize. [number_chains(i,j) = 0] is used as a flag below.
if ( method == CONTINUATION_METHOD ) {
    number_chains.resize(rows,cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            number_chains(i,j) = 0;
        }
    }

    chains.resize(rows,cols);
    chain_edges.resize(rows,cols);
//    iplus.resize(rows,cols);
//    iminus.resize(rows,cols);
//    jplus.resize(rows,cols);
//    jminus.resize(rows,cols);
    chain_list.clear();
}

    std::vector<int> coordinates(2);

    /* Recall that the order of vertex index on the rectangles is given by hcube;
       it is not trivial, the order is:
                       3-vertex +---+ 2-vertex
                                |\  |
                                | \ |
                                |  \|
                       1-vertex +---+ 0-vertex

       The order of the faces is also given by hcube:
                                4-face
                                +----+
                                |\   |
                        2-face  | \  |  3-face
                                |  \ |
                                +----+
                                0-face
       The diagonal connectic 0-vertex to 3-vertex is 1-face.
    */

//    std::ofstream off("cell.txt", std::ifstream::out);

    // loop over the rectangles
    for (int i = 0; i < rows; i++) { //std::cout << "i = " << i << std::endl;
        for (int j = 0; j < cols; j++) { //std::cout << "    j = " << j << std::endl;
            vert[0*hn + 0] = gv->grid(i + 1, j + 0).component(0); // u + du;
            vert[1*hn + 0] = gv->grid(i + 0, j + 0).component(0); // u;
            vert[2*hn + 0] = gv->grid(i + 1, j + 1).component(0); // u + du;
            vert[3*hn + 0] = gv->grid(i + 0, j + 1).component(0); // u;

            vert[0*hn + 1] = gv->grid(i + 1, j + 0).component(1); // v;
            vert[1*hn + 1] = gv->grid(i + 0, j + 0).component(1); // v;
            vert[2*hn + 1] = gv->grid(i + 1, j + 1).component(1); // v + dv;
            vert[3*hn + 1] = gv->grid(i + 0, j + 1).component(1); // v + dv;

//            std::cout << "    After verts" << std::endl;

//            off << "    Cell " << i << ", " << j << " is: " << gv->cell_type(i, j) << std::endl;

            if (gv->cell_type(i, j) != CELL_IS_INVALID) {
                    if (gv->cell_type(i, j) == CELL_IS_SQUARE){
                        ssimp = nsimp_;
                        sface = nface_;
                    } 
                    else if (gv->cell_type(i, j) == CELL_IS_TRIANGLE) {
                        ssimp = tsimp;
                        sface = tface;
                    }
                    else continue;

                    // The posibility of a zero on a cell is given by this...
                    zero = 0; // isso vai substituir o usso de uma variavel logica (false) TODO : --> colocar como booleano (zero --> root)

                    // foncub is filled here
                    int status = impf->function_on_square(&foncub[0][0], i, j);
                    
                    if (status != 0) {
                        refval = foncub[0][0];
                        for (int l = 1; l < 4; l++) {
                            if (l == 2 && gv->cell_type(i, j) == CELL_IS_TRIANGLE) {
                            }
                            else {
                                if (refval * foncub[0][l] < 0.0) zero = 1;
                            }
                        }
                    }

                    //std::cout << "    zero = " << zero << std::endl;

                    if (zero != 0) {
                        //  solve for the current cube
                        double u[hn][hm + 1]; //TODO Remove
                        double g[hm][hm+1];
                        double stormd[hm];

                        nsoln_ = hc.cubsol(solptr_, sol_, dims_, sptr_, nsoln_, &foncub[0][0],
                                           exstfc, face_, facptr_, dimf_, vert,
                                           ncvert_, hn, hm, ssimp, nsface_, sface, &u[0][0], &g[0][0],
                                           stormd, storm_);

                        // IMPROVE THE SOLUTION USING A ZERO-FINDER.
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

                                sol_[0*dims_ + sp] = p_improved_newton.component(0); // sol(1,sp) = v(1)
                                sol_[1*dims_ + sp] = p_improved_newton.component(1); // sol(2,sp) = v(2)

                            }
                        }

                        // MAKE THE LIST OF EDGE POINTERS
                        int N_EDGES = hc.mkedge(edges_, dime_, nedges_, smpedg_, solptr_,
                                            fnbr_, ssimp, nsface_); // TODO: Em Fortran esta nsimp_

                        nedges_ = N_EDGES;

                        //STORE ALL PAIRS OF EDGES THAT WERE FOUND
                        if (nedges_ > 0) {
                            // A Topological Sort is called for ordering the segments in each cell.
                            //
                            if ( method == CONTINUATION_METHOD ) {
                                if( !topological_sort(i,j) ) continue;
                                coordinates[0] = i; coordinates[1] = j;
                                // Store the coodinates of the GridValues that have solution, so after
                                // the sort, a list of order n is used and not the whole grid n^2.
                                chain_list.push_back(coordinates);
//                                cout << "Coordinates: " << coordinates[0] << ", " << coordinates[1] << endl;
//                                cout << "Number of chains: " << number_chains(i,j) << endl;
//                                cout << "The first point:  " << chains(i,j)[0][0] << endl;
                            }

                            // The segments are store in pairs
                            //
                            if ( method == SEGMENTATION_METHOD ) {
                                for (nedg = 0; nedg < nedges_; nedg++) {

                                   // Store the segments
                                   RealVector p1(2), p2(2);
                                   p1.component(0) = sol_[0*dims_ + edges_[0*dime_ + nedg ]];
                                   p1.component(1) = sol_[1*dims_ + edges_[0*dime_ + nedg ]];

                                   // TODO: Pablo comentou a seguinte linha em 11 Janeiro 2012
                                   //       if (edges_[1*dime_ + nedg ] < 0) return 0;
                                   p2.component(0) = sol_[0*dims_ + edges_[1*dime_ + nedg ]];
                                   p2.component(1) = sol_[1*dims_ + edges_[1*dime_ + nedg ]];

                                   vrs.push_back(p1);
                                   vrs.push_back(p2);
                                }
                            }
                        }
                    }

                    // 
            }
        }
    }

//    off.close();

    if ( method == SEGMENTATION_METHOD ) {
//        printf("ENDS:     Contour2D()/For Segments\n\n");
        return 0;
    }

// TEST
//   cout << "Cells: ";
//   for (int i = 0; i < chain_list.size(); i++) {
//       cout << number_chains(chain_list[i][0],chain_list[i][1]) << " ";
//   }
//   cout << endl;
// End TEST

    // The following integers keep the path of the concatenation
    int prev, middle, next;

    curves.clear();
    is_circular.clear();

    // The loop is given over the active cell grid stored in [chain_list]
    for (int i = 0; i < chain_list.size(); /*For two segments, the cell is used twice*/) {
        int i_index = chain_list[i][0]; int i_start = i_index;
        int j_index = chain_list[i][1]; int j_start = j_index;

        // Here is the increment of [int i], when needed
        if ( number_chains(i_index, j_index) < 2 ) i++;
        if ( number_chains(i_index, j_index) == 0 ) continue;

        // Reused for each branch, is a <deque> because it can be filled by both sides.
        std::deque <RealVector> branch;
        bool circular_is = false;
        int chain_size = chains(i_index,j_index)[0].size();

        // [prev] store the direction not followed, [middle] the direction to follow; [next] will
        // store the entering direction in the next cell.
        prev   = chain_edges(i_index,j_index)[0][0];
        middle = chain_edges(i_index,j_index)[0][1];

        // Store the points in the [branch] and eliminate the points of the list 
        for (int j = 0; j < chain_size; j++) branch.push_back( chains(i_index,j_index)[0][j] );
        chains(i_index,j_index)[0].erase( chains(i_index,j_index)[0].begin());
        chain_edges(i_index,j_index)[0].erase( chain_edges(i_index,j_index)[0].begin());

        do {
            restore:

            // [middle] also works as a flag. If it is used once then it is setted as -1, there are
            //          no previous direction, then returning to the same cell end the procedure.
            if (middle == -1) break;
            number_chains(i_index, j_index)--;

            // [middle] gives the direction to follow and the index of the next point
            switch ( middle ) {
                case 0: {
                    j_index--;
                    next = 4;
                }
                break;

                case 2: {
                    i_index--;
                    next = 3;
                }
                break;

                case 3: {
                    i_index++;
                    next = 2;
                }
                break;

                case 4: {
                    j_index++;
                    next = 0;
                }
                break;

                default:
                break;
            }
            // Preventing to leave the domain.
            if ( (i_index < 0) || (j_index < 0) || (i_index >= rows) || (j_index >= cols) || middle == 1 ) {
                middle = prev;
                i_index = i_start;
                j_index = j_start;
                prev = -1;
//                /* DEBUG... */
//                cout << " :: FIZ __RESTORE__" << endl;
//                cout << "index(" << next << "): " << i_index << ", " << j_index << endl;
//                cout << "i: " << i_index << ", " << i_start << endl;
//                cout << "j: " << j_index << ", " << j_start << endl;
//                /* ... END DEBUG */
                goto restore;
            }

            // Some restriction may cut the continuation of the curve, the following IF saves this issues.
            if( chain_edges(i_index,j_index).size() == 0 ) break;
            // Store the points in the [branch] and eliminate the points of the list.
            // The chain in (i, j) must be identified with [next] the direction of push_back or
            // push_front is denoted by the use of [prev = -1], else it will be the middle direction.
            for (int k = 0; k < chains(i_index,j_index)[0].size(); k++) {
                if ( chain_edges(i_index,j_index)[k][0] == next ) {
                    chain_size = chains(i_index,j_index)[k].size();
                    if (prev == -1)
                        for (int j = 1; j < chain_size; j++) branch.push_front( chains(i_index,j_index)[k][j] );
                    else
                        for (int j = 1; j < chain_size; j++) branch.push_back( chains(i_index,j_index)[k][j] );
                    middle = chain_edges(i_index,j_index)[k][1];
                    chains(i_index,j_index).erase( chains(i_index,j_index).begin() + k );
                    chain_edges(i_index,j_index).erase( chain_edges(i_index,j_index).begin() + k);
                    if ( (i_index == i_start) && (j_index == j_start) && (next == prev) ) {
                        // TODO: The repeated point it is not deleted! It is convenient to keep it.
                        // branch.erase( branch.begin() );
                        circular_is = true;
                    }
                    break;
                }
                if ( chain_edges(i_index,j_index)[k][1] == next ) {
                    chain_size = chains(i_index,j_index)[k].size();
                    if (prev == -1) 
                        for (int j = chain_size-2; j >= 0; j--) branch.push_front( chains(i_index,j_index)[k][j] );
                    else            
                        for (int j = chain_size-2; j >= 0; j--) branch.push_back( chains(i_index,j_index)[k][j] );
                    middle = chain_edges(i_index,j_index)[k][0];
                    chains(i_index,j_index).erase( chains(i_index,j_index).begin() + k );
                    chain_edges(i_index,j_index).erase( chain_edges(i_index,j_index).begin() + k);
                    if ( (i_index == i_start) && (j_index == j_start) && (next == prev) ) {
                        // TODO: The repeated point it is not deleted! It is convenient to keep it.
                        // branch.erase( branch.begin() );
                        circular_is = true;
                    }
                    break;
                }
            }

        } while ( number_chains(i_index, j_index) > 0 );

        curves.push_back(branch);
        is_circular.push_back(circular_is);
    }

// For DEBUGING, the following lines can be included instead of the two ereasing statements.
///* DEBUG...
//cout << "Antes de deletar tenho: "<< endl;
//for(int ii = 0; ii < chains(i_index,j_index).size(); ii++){
//    for (int jj = 0; jj < chains(i_index,j_index)[ii].size(); jj++){
//        cout << " " << chains(i_index,j_index)[ii][jj];
//    }
//    cout << endl;
//}
//... END DEBUG */
//                    chains(i_index,j_index).erase( chains(i_index,j_index).begin() + k );
//                    chain_edges(i_index,j_index).erase( chain_edges(i_index,j_index).begin() + k);
///* DEBUG...
//cout << "Depois de deletar tenho: "<< endl;
//for(int ii = 0; ii < chains(i_index,j_index).size(); ii++){
//    for (int jj = 0; jj < chains(i_index,j_index)[ii].size(); jj++){
//        cout << " " << chains(i_index,j_index)[ii][jj];
//    }
//    cout << endl;
//}
//... END DEBUG */

///*DEBUG: Impressions*/
//    cout << "Finalizei com " << curves.size() << " curvas :: ::   ";
//    for (int k = 0; k < curves.size(); k++) cout << "Tamanho de " << k << " es: " << curves[k].size() << " y es circular: " << is_circular[k] << "; ";
//    cout << endl;

//    for (int k = 0; k < curves.size(); k++) {
//        cout << "Curva " << k << ":";
//        for (int l = 0; l < curves[k].size(); l++) cout << " " << curves[k][l];
//        cout << endl;
//    }
///*END DEBUG...*/

    printf("ENDS:     Contour2D()/Continuous\n\n");

    return 0;
}


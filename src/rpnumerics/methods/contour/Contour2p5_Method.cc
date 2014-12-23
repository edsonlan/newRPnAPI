#include "Contour2p5_Method.h"

// Static variables are defined here:
//
bool Contour2p5_Method::is_first = true;

HyperCube Contour2p5_Method::hc;

int Contour2p5_Method::hn = 3;
int Contour2p5_Method::hm = 2;

// Estos numeros todavia no fueron determindados...
int Contour2p5_Method::nsface_;  // =  3; colocarlos en la primera salida de impresion y ver cuanto valen.
int Contour2p5_Method::nface_;   // =  5; TODO: pero no es importante.
int Contour2p5_Method::nsoln_  = -1;
int Contour2p5_Method::nedges_;
// Hasta aqui. (Los dos ultimos no parecen ser importantes.)

int Contour2p5_Method::dims_   = 18;
int Contour2p5_Method::dime_   = 36;
int Contour2p5_Method::dimf_   = 18;
int Contour2p5_Method::ncvert_ =  8;
int Contour2p5_Method::dncv    =  8; //TODO, ver si este no es el ncvert_
int Contour2p5_Method::nsimp_  =  6;

int Contour2p5_Method::numberOfCombinations = 4;

int * Contour2p5_Method::storn_;
int * Contour2p5_Method::storm_;

Matrix<double> Contour2p5_Method::cvert_;
Matrix<double> Contour2p5_Method::vert;
Matrix<int>    Contour2p5_Method::bsvert_;
Matrix<int>    Contour2p5_Method::perm_;
Matrix<int>    Contour2p5_Method::comb_;
Matrix<double> Contour2p5_Method::foncub;
Matrix<int>    Contour2p5_Method::fnbr_;
Matrix<int>    Contour2p5_Method::facptr_;
Matrix<int>    Contour2p5_Method::face_;
Matrix<double> Contour2p5_Method::cpp_sol;
Matrix<int>    Contour2p5_Method::solptr_;
Matrix<int>    Contour2p5_Method::cpp_edges_;
Matrix<int>    Contour2p5_Method::smpedg_;
int * Contour2p5_Method::exstfc;
int * Contour2p5_Method::sptr_;

//int Contour2p5_Method::tsimp = 1;
//int Contour2p5_Method::tface = 3;

int * Contour2p5_Method::index;
int * Contour2p5_Method::usevrt;

void Contour2p5_Method::allocate_arrays(void){
    if (is_first){
        storn_ = new int[hn + 1];
        storm_ = new int[hm + 1];
                
        cvert_.resize(ncvert_, hn); // Was transposed: hccube.inc:4
        vert.resize(ncvert_, hn);
        bsvert_.resize(hn + 1, hn); // Was transposed: hccube.inc:5
        perm_.resize(hn, nsimp_); // Was NOT transposed: hccube.inc:5
        comb_.resize(numberOfCombinations, hm + 1); // Was transposed: hccube.inc:7
        foncub.resize(hm, ncvert_); // NOT transposed, as seen in hcsoln.inc:3.
        nsface_ = hc.mkcomb(comb_.data(), hn + 1, hm + 1);
        fnbr_.resize(nsface_, nsface_);
        cpp_sol.resize(hn, dims_); // NOT transposed, as seen in hcsoln.inc:3.
        solptr_.resize(nsimp_, nsface_); // Transposed, as seen in hcsoln.inc:4.
        cpp_edges_.resize(2, dime_); // NOT transposed as seen in hcedge.inc:3.
        smpedg_.resize(nsimp_, 2); // Transposed, as seen in hcedge.inc:3.
        for (int i = 0; i < nsimp_*2; i++) smpedg_(i) = 0; //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos
        facptr_.resize(nsimp_, nsface_);
        face_.resize(hm + 1, dimf_); // Was NOT transposed: hccube.inc:7
                                    // Confirmed in: hcmarc.inc:4.

        hc.mkcube(cvert_.data(), bsvert_.data(), perm_.data(), ncvert_, nsimp_, hn);

        nface_ = hc.mkface(face_.data(), facptr_.data(), fnbr_.data(), dimf_, nsimp_, hn, hm, nsface_,
                           bsvert_.data(), comb_.data(), perm_.data(), &storn_[0], &storm_[0]);

        exstfc = new int[nface_]; for (int i = 0; i < nface_; i++) exstfc[i] = 1; // Verify that nface_ is DNFACE.
        sptr_  = new int[nface_];

        index = new int[4]; index[0] = 0; index[1] = 2; index[2] = 3; index[3] = 1; 
        
        usevrt = new int[dncv]; // Defined in: contour3.F.
        for (int i = 0; i < dncv; i++) usevrt[i] = 1;

        is_first = false;

        printf("++++++++++++++++ Contour2p5_Method: REMEMBER TO INVOKE deallocate_arrays() AT QUIT-TIME!!!\n++++++++++++++++ DON\'T SAY I DIDN\'T WARN YOU!!!\n");
        printf("    After allocating arrays, nsface_ = %d, nface_ =  %d\n", nsface_, nface_);
    }

    return;
}

void Contour2p5_Method::deallocate_arrays(void){
    if (!is_first){
        delete [] usevrt;
        delete [] index;
    
        delete [] sptr_;
        delete [] exstfc;

        delete [] storm_;
        delete [] storn_;

        is_first = true;
    }

    printf("++++++++++++++++ Contour2p5_Method: arrays deallocated. ++++++++++++++++\n");

    return;
}

/* HERE */
void Contour2p5_Method::contour2p5(TwoImplicitFunctions *timpf, std::vector<RealVector> &curve_vrs, std::vector<RealVector> &domain_vrs){
    allocate_arrays();

    curve_vrs.clear();
    domain_vrs.clear(); 
    
    cout<<"Tamanho de oc : "<<timpf->segment_value()->size()<<endl;

    for (int i = 0; i < timpf->segment_value()->size()/2; i++) {
        if ( !timpf->valid_segment(2*i) ){
            cout<<"Segmento "<<i<<" invalido"<<endl;
            continue;
        }
            
            

        curve2p5(timpf, 2*i, curve_vrs, domain_vrs);
    }
    
    return;
}

void Contour2p5_Method::contour2p5_for_curve(TwoImplicitFunctions *timpf, std::vector<RealVector> &curve_vrs, std::vector<RealVector> &domain_vrs){
    allocate_arrays();

    curve_vrs.clear();
    domain_vrs.clear();

    for (int i = 0; i < timpf->segment_value()->size() - 1; i++) {
        if ( !timpf->valid_segment(i) ) continue;
        curve2p5(timpf, i, curve_vrs, domain_vrs);
    }
    
    return;
}

// TODO: Each point is being extended twice: it belongs to two segments.
//       This situation should be remedied, by extending not a segment, but a point individually.
//
void Contour2p5_Method::curve2p5(TwoImplicitFunctions *timpf, int current_segment_index, 
                                 std::vector<RealVector> &curve_vrs, 
                                 std::vector<RealVector> &domain_vrs){
                                 
    // Get the current segment, etc.
    //
    GridValues *gv = timpf->grid_value();
    std::vector<RealVector> *segments = timpf->segment_value();
                                                                  
    std::vector<RealVector> current_segment;
    current_segment.resize(2);
    for (int j = 0; j < 2; j++){
        //current_segment[j].resize(2);
        //for (int k = 0; k < 2; k++) current_segment[j].component(k) = segments->at(current_segment_index + j).component(k);
        current_segment[j] = segments->at(current_segment_index + j);
    }

    // Get the current coordinates to be used later.
    //                                 
    double ub = current_segment[0].component(0);
    double vb = current_segment[0].component(1);
    double dub = current_segment[1].component(0) - ub;
    double dvb = current_segment[1].component(1) - vb;

    double work1[hn][hm + 1]; //  Equivalent to double u[hn][hm + 1]; Was NOT transposed, as seen in hccube.inc:10.
    double work2[hm][hm + 1]; //  Equivalent to double g[hm][hm + 1]; Was 
    double stormd[hm];        //  Equivalent to double stormd[hm];

    for (int i = 0; i < gv->grid.rows() - 1; i++) {
        for (int j = 0; j < gv->grid.cols() - 1; j++) {
            if (gv->cell_type(i, j) != CELL_IS_SQUARE) continue; // Only for squares within the domain. Triangles are not dealt with yet.

            // Check for duplicate squares in singular case
            if ( timpf->is_singular() ){
                //c to do: import the left sizes for testing closeness
                double dumax = gv->grid_resolution.component(0)*2.0;
                double dvmax = gv->grid_resolution.component(1)*2.0;

                double bmid[2];
                bmid[0] = ub + .5 * dub;
                bmid[1] = vb + .5 * dvb;

                double dur = gv->grid_resolution[0];
                double dvr = gv->grid_resolution[1];

                // Equivalent to Double_Contact::left_right_adjacency().
                //            grid(i, j) is used to set the middle cell point.
                if ((fabs(gv->grid(i, j).component(0) + .5 * dur - bmid[0]) <= dumax) &&
                    (fabs(gv->grid(i, j).component(1) + .5 * dvr - bmid[1]) <= dvmax)) continue;
            }

            // function_on_cube will call the explicit function_on_vertices for each method
            //
            int info = function_on_cube(timpf, foncub, i, j); 
            if (info == VALID_FUNCTION_ON_VERTICES) {
                nsoln_ = hc.cpp_cubsol(solptr_.data(), cpp_sol, dims_,
                        &sptr_[0], nsoln_, foncub.data(), &exstfc[0],
                        face_.data(), facptr_.data(), dimf_, cvert_.data(),
                        ncvert_, hn, hm, nsimp_, nsface_, nface_, &work1[0][0],
                        &work2[0][0], &stormd[0], &storm_[0]);

                // Make the list of edge pointers.
                nedges_ = hc.cpp_mkedge(cpp_edges_, dime_, nedges_, smpedg_.data(),
                        solptr_.data(), fnbr_.data(), nsimp_, nsface_);

                // Fill the segments both on domain and on curve
                filedg3(cpp_sol, cpp_edges_,
                    nedges_, gv->grid_resolution,
                    current_segment,
                    gv->grid(i, j),
                    curve_vrs, domain_vrs);
            }
        } // For j
    } // For i

    return;
}

void Contour2p5_Method::filedg3(Matrix<double> &sol_, Matrix<int> &edges_,
                                int nedges_,
                                const RealVector &grid_resolution,
                                const std::vector<RealVector> &segment,
                                const RealVector &uv,
                                std::vector<RealVector> &curve_segments, std::vector<RealVector> &domain_segments) {

    // Set base points, increments, and midpoint
    //
    double ub = segment[0].component(0);
    double vb = segment[0].component(1);
    double dub = segment[1].component(0) - ub;
    double dvb = segment[1].component(1) - vb;

    double u = uv.component(0); // ur0 + dur * (ir-1)
    double v = uv.component(1); // vr0 + dvr * (jr-1)

    double dur = grid_resolution.component(0);
    double dvr = grid_resolution.component(1);

    // Store all pairs of edges that were found
    RealVector p1(2), p2(2);
    for (int nedg = 0; nedg < nedges_; nedg++) {
        // For the curve
        p1.component(0) = ub + dub * sol_(2, edges_(0, nedg)); // LX1 Was  = segend(1,1,sn) = ub + dub * sol(3,edges(1,nedg))
        p1.component(1) = vb + dvb * sol_(2, edges_(0, nedg)); // LY1 Was  = segend(2,1,sn) = vb + dvb * sol(3,edges(1,nedg))

        p2.component(0) = ub + dub * sol_(2, edges_(1, nedg)); // LX2 Was  = segend(1,2,sn) = ub + dub * sol(3,edges(2,nedg))
        p2.component(1) = vb + dvb * sol_(2, edges_(1, nedg)); // LY2 Was  = segend(2,2,sn) = vb + dvb * sol(3,edges(2,nedg))

        curve_segments.push_back(p1);
        curve_segments.push_back(p2);

        // For the domain;
        p1.component(0) = u + dur * sol_(0, edges_(0, nedg)); // RX1 Was:  = segend(1,1,sn) = u + dur * sol(1,edges(1,nedg))
        p1.component(1) = v + dvr * sol_(1, edges_(0, nedg)); // RY1 Was:  = segend(2,1,sn) = v + dvr * sol(2,edges(1,nedg))

        p2.component(0) = u + dur * sol_(0, edges_(1, nedg)); // RX2 Was:  = segend(1,2,sn) = u + dur * sol(1,edges(2,nedg))
        p2.component(1) = v + dvr * sol_(1, edges_(1, nedg)); // RY2 Was:  = segend(2,2,sn) = v + dvr * sol(2,edges(2,nedg))

        domain_segments.push_back(p1);
        domain_segments.push_back(p2);
    }

    return;
}


int Contour2p5_Method::function_on_cube(TwoImplicitFunctions *timpf, Matrix<double> &foncub, int ir, int jr) {
    bool zero[2] = {false, false};

    double val[2];
    double refval[2];

    if (timpf->function_on_vertices(refval, ir, jr, 0) == INVALID_FUNCTION_ON_VERTICES) return INVALID_FUNCTION_ON_VERTICES;

    for (int kl = 0; kl < 2; kl++) {
        for (int kr = 0; kr < 4; kr++) {
            int domain_i, domain_j;

            // The domain index is given by ordenation inside HyperCube.
            /**
             *    kr = 3  +-----+  kr = 2
             *            + \   +
             *            +  \  +
             *            +   \ +
             *    kr = 0  +-----+  kr = 1
            **/
            if (kr == 0) {
                domain_i = ir;
                domain_j = jr;
            } else if (kr == 1) {
                domain_i = ir + 1;
                domain_j = jr;
            } else if (kr == 2) {
                domain_i = ir + 1;
                domain_j = jr + 1;
            } else if (kr == 3) {
                domain_i = ir;
                domain_j = jr + 1;
            }

            if (timpf->function_on_vertices(val, domain_i, domain_j, kl) == INVALID_FUNCTION_ON_VERTICES) return INVALID_FUNCTION_ON_VERTICES;

            for (int comp = 0; comp < 2; comp++) {
                foncub(comp, kl + 2 * index[kr]) = val[comp];

                if (refval[comp] * val[comp] < 0.0) zero[comp] = true;
            }
        }
    }

    if (!zero[0] && !zero[1]) return NO_ZERO_ON_CUBE;

    return VALID_FUNCTION_ON_VERTICES;

}


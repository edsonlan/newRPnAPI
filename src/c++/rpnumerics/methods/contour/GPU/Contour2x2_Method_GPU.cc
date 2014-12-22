#include "Contour2x2_Method_GPU.h"

// Static variables are defined here:
//
bool Contour2x2_Method_GPU::is_first = true;

HyperCube Contour2x2_Method_GPU::hc;

int Contour2x2_Method_GPU::hn = 4;
int Contour2x2_Method_GPU::hm = 3;

// Estos numeros todavia no fueron determindados...
int Contour2x2_Method_GPU::nsface_;  // =  3; colocarlos en la primera salida de impresion y ver cuanto valen.
int Contour2x2_Method_GPU::nface_;   // =  5;TODO: pero no es importante.
int Contour2x2_Method_GPU::nsoln_  = -1;
int Contour2x2_Method_GPU::nedges_;
// Hasta aqui. (Los dos ultimos no parecen ser importantes.)

int Contour2x2_Method_GPU::dims_   = 125;
int Contour2x2_Method_GPU::dime_   = 200;
int Contour2x2_Method_GPU::dimf_   =  84;
int Contour2x2_Method_GPU::ncvert_ =  16;
int Contour2x2_Method_GPU::dncv    =  16; //TODO, ver si este no es el ncvert_
int Contour2x2_Method_GPU::nsimp_  =  24;

int Contour2x2_Method_GPU::numberOfCombinations = 5;

double Contour2x2_Method_GPU::ul0;
double Contour2x2_Method_GPU::vl0;
double Contour2x2_Method_GPU::ur0;
double Contour2x2_Method_GPU::vr0;

double Contour2x2_Method_GPU::dul;
double Contour2x2_Method_GPU::dvl;
double Contour2x2_Method_GPU::dur;
double Contour2x2_Method_GPU::dvr;

double Contour2x2_Method_GPU::dumax;
double Contour2x2_Method_GPU::dvmax;

int * Contour2x2_Method_GPU::storn_;
int * Contour2x2_Method_GPU::storm_;

Matrix<double> Contour2x2_Method_GPU::cvert_;
Matrix<double> Contour2x2_Method_GPU::vert;
Matrix<int>    Contour2x2_Method_GPU::bsvert_;
Matrix<int>    Contour2x2_Method_GPU::perm_;
Matrix<int>    Contour2x2_Method_GPU::comb_;
Matrix<double> Contour2x2_Method_GPU::foncub;
Matrix<int>    Contour2x2_Method_GPU::fnbr_;
Matrix<int>    Contour2x2_Method_GPU::facptr_;
Matrix<int>    Contour2x2_Method_GPU::face_;
Matrix<double> Contour2x2_Method_GPU::cpp_sol;
Matrix<int>    Contour2x2_Method_GPU::solptr_;
Matrix<int>    Contour2x2_Method_GPU::cpp_edges_;
Matrix<int>    Contour2x2_Method_GPU::smpedg_;
int * Contour2x2_Method_GPU::exstfc;
int * Contour2x2_Method_GPU::sptr_;

//int Contour2x2_Method_GPU::tsimp = 1;
//int Contour2x2_Method_GPU::tface = 3;

int * Contour2x2_Method_GPU::index;
int * Contour2x2_Method_GPU::usevrt;

void Contour2x2_Method_GPU::allocate_arrays(void){
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
        //TODO: Ver o que acontece, pois se nao sao inicializadas coloca valores estranhos
        for (int i = 0; i < nsimp_*2; i++) smpedg_(i) = 0;
        facptr_.resize(nsimp_, nsface_);
        face_.resize(hm + 1, dimf_); // Was NOT transposed: hccube.inc:7
                                     // Confirmed in: hcmarc.inc:4.

        hc.mkcube(cvert_.data(), bsvert_.data(), perm_.data(), ncvert_, nsimp_, hn);

//// DEBUG perm e bsvert
//    cout << "perm(" << hn << ", " << nsimp_ << "):" << endl;
//    for (int i = 0; i < hn; i++) {
//        for (int j = 0; j < nsimp_; j++) {
//            cout << " " << perm_(i,j);
//        }
//        cout << endl;
//    }
//    cout << "bsvert(" << hn << ", " << hn+1 << "):" << endl;
//    for (int i = 0; i < hn; i++) {
//        for (int j = 0; j < hn+1; j++) {
//            cout << " " << bsvert_(j,i);
//        }
//        cout << endl;
//    }
//// END DEBUG

        nface_ = hc.mkface(face_.data(), facptr_.data(), fnbr_.data(), dimf_, nsimp_, hn, hm, nsface_,
                           bsvert_.data(), comb_.data(), perm_.data(), &storn_[0], &storm_[0]);

        // Verify that nface_ is DNFACE.
        exstfc = new int[nface_]; for (int i = 0; i < nface_; i++) exstfc[i] = 1;
        sptr_  = new int[nface_];

        index = new int[4]; index[0] = 0; index[1] = 2; index[2] = 3; index[3] = 1; 
        
        usevrt = new int[dncv]; // Defined in: contour3.F.
        for (int i = 0; i < dncv; i++) usevrt[i] = 1;

        is_first = false;

        printf("++++++++++++++++ Contour2x2_Method_GPU: REMEMBER TO INVOKE deallocate_arrays() AT QUIT-TIME!!!\n++++++++++++++++ DON\'T SAY I DIDN\'T WARN YOU!!!\n");
        printf("    After allocating arrays, nsface_ = %d, nface_ =  %d\n", nsface_, nface_);
    }

    return;
}

void Contour2x2_Method_GPU::deallocate_arrays(void){
    if (!is_first){
        delete [] usevrt;
        delete [] index;
    
        delete [] sptr_;
        delete [] exstfc;

        delete [] storm_;
        delete [] storn_;

        is_first = true;
    }

    printf("++++++++++++++++ Contour2x2_Method_GPU: arrays deallocated. ++++++++++++++++\n");

    return;
}

void Contour2x2_Method_GPU::curve2x2(ThreeImplicitFunctions *timpf,
                                 std::vector<RealVector> &left_vrs,   // on_domain
                                 std::vector<RealVector> &right_vrs){ // on_curve
                            
    allocate_arrays();
     
    // Get the current data.
    //
    GridValues *gv_left  = timpf->grid_value_left();
    GridValues *gv_right = timpf->grid_value_right();

    ul0 = gv_left->grid(0,0).component(0);
    vl0 = gv_left->grid(0,0).component(1);
    ur0 = gv_right->grid(0,0).component(0);
    vr0 = gv_right->grid(0,0).component(1);

    dul = gv_left->grid_resolution.component(0);
    dvl = gv_left->grid_resolution.component(1);
    dur = gv_right->grid_resolution.component(0);
    dvr = gv_right->grid_resolution.component(1);

    dumax = 2.0 * max ( dur, dul );
    dvmax = 2.0 * max ( dvr, dvl );

    double u[hn][hm + 1];
    double g[hm][hm + 1];
    double stormd[hm];

    for (int il = 0; il < gv_left->grid.rows() - 1; il++) {
        for (int jl = 0; jl < gv_left->grid.cols() - 1; jl++) {

            // Only for squares within the domain.
            if (gv_left->cell_type(il, jl) != CELL_IS_SQUARE) {
                continue;
            
            }
//            if ( !(gv_left->cell_is_real(il, jl)) ) continue;


            if(!timpf->prepare_cell(il, jl)) continue;

            for (int ir = 0; ir < gv_right->grid.rows() - 1; ir++){
                for (int jr = 0; jr < gv_right->grid.cols() - 1; jr++){
                    if (gv_right->cell_type(ir, jr) == CELL_IS_SQUARE){
                        if ( (timpf->is_singular()) && left_right_adjacency(il, jl, ir, jr)) continue;

                        if (filhcub4(timpf, ir, jr, index, foncub.data())) {

                            nsoln_ = hc.cpp_cubsol(solptr_.data(), cpp_sol, dims_, 
                                                   &sptr_[0], nsoln_, foncub.data(), &exstfc[0], 
                                                   face_.data(), facptr_.data(), dimf_, cvert_.data(), 
                                                   ncvert_, hn, hm, nsimp_, nsface_, nface_, &u[0][0], 
                                                   &g[0][0], &stormd[0], &storm_[0]);
                        
                            nedges_ = hc.cpp_mkedge(cpp_edges_, dime_, nedges_, smpedg_.data(), 
                                                    solptr_.data(), fnbr_.data(), nsimp_, nsface_);
                           
                            filedg4(cpp_sol, cpp_edges_, nedges_,
                                     il, jl, ir, jr, left_vrs, right_vrs);

                        }
                    }
                } // For jr
            } // For ir
        } // For jl
    } // For il

    return;
}

void Contour2x2_Method_GPU::cpu_curve2x2(ThreeImplicitFunctions *timpf,
                                     const std::vector<short int> &cells_left,
                                     const std::vector<short int> &cells_right,
                                     std::vector<RealVector> &left_vrs,   // on_domain
                                     std::vector<RealVector> &right_vrs){ // on_curve

    allocate_arrays();
     
    // Get the current data.
    //
    GridValues *gv_left  = timpf->grid_value_left();
    GridValues *gv_right = timpf->grid_value_right();

    ul0 = gv_left->grid(0,0).component(0);
    vl0 = gv_left->grid(0,0).component(1);
    ur0 = gv_right->grid(0,0).component(0);
    vr0 = gv_right->grid(0,0).component(1);

    dul = gv_left->grid_resolution.component(0);
    dvl = gv_left->grid_resolution.component(1);
    dur = gv_right->grid_resolution.component(0);
    dvr = gv_right->grid_resolution.component(1);

    dumax = 2.0 * max ( dur, dul );
    dvmax = 2.0 * max ( dvr, dvl );

    double u[hn][hm + 1];
    double g[hm][hm + 1];
    double stormd[hm];

    for (int running_index = 0; running_index < cells_left.size(); running_index++) {
        int il = std::floor(cells_left[running_index]/gv_left->grid.cols());
        int jl = cells_left[running_index] - il*gv_left->grid.cols();

        // Only for squares within the domain.
        // TODO: This must change in order to deal with the hypothenuse of the boundary.
        //
        if (gv_left->cell_type(il, jl) != CELL_IS_SQUARE) {
            continue;
        }

        if(!timpf->prepare_cell(il, jl)) continue;

        int ir = std::floor(cells_right[running_index]/gv_right->grid.cols());
        int jr = cells_right[running_index] - ir*gv_right->grid.cols();
        
        if (gv_right->cell_type(ir, jr) == CELL_IS_SQUARE){
            if ( (timpf->is_singular()) && left_right_adjacency(il, jl, ir, jr)) continue;

            if (filhcub4(timpf, ir, jr, index, foncub.data())) {

                nsoln_ = hc.cpp_cubsol(solptr_.data(), cpp_sol, dims_, 
                                                   &sptr_[0], nsoln_, foncub.data(), &exstfc[0], 
                                                   face_.data(), facptr_.data(), dimf_, cvert_.data(), 
                                                   ncvert_, hn, hm, nsimp_, nsface_, nface_, &u[0][0], 
                                                   &g[0][0], &stormd[0], &storm_[0]);
                        
                nedges_ = hc.cpp_mkedge(cpp_edges_, dime_, nedges_, smpedg_.data(), 
                                                    solptr_.data(), fnbr_.data(), nsimp_, nsface_);
                           
                filedg4(cpp_sol, cpp_edges_, nedges_,
                                 il, jl, ir, jr, left_vrs, right_vrs);

            }
        }
    }

    return;
}

bool Contour2x2_Method_GPU::filhcub4(ThreeImplicitFunctions *timpf,
                                 int ir, int jr, int *index, double *foncub){
    bool zero[3] = {false, false, false};

    double val[3];    // To be filled e.g. by Double_Contact::function_on_cell();
    double refval[3]; // To be filled e.g. by Double_Contact::function_on_cell();
    
    if (!timpf->function_on_cell(refval, ir, jr, 0, 0)) return false;

    for (int kl = 0; kl < 4; kl++){
        for (int kr = 0; kr < 4; kr++){
            if (!timpf->function_on_cell(val, ir, jr, kl, kr)) return false;

            for (int comp = 0; comp < 3; comp++){
                foncub[comp*ncvert_ + 4*index[kl] + index[kr]] = val[comp];
                // Modified by Morante on 21-06-2011 by advice from Castaneda.
                if (refval[comp]*val[comp] < 0.0) zero[comp] = true;
            }
        }
    }
//    // DEBUG: Sufficient condition:
//    //        Probably next line must be increased becaus index[2] = 3:
//    //        if (!timpf->function_on_cell(val, ir, jr, 2, 2)) return false;
//    if ( (refval[0]*val[0] < 0.0) && (refval[1]*val[1] < 0.0) && (refval[2]*val[2] < 0.0) ) {
//        cout << endl;
//        cout << "***** Sufficient (" << ir << ", " << jr << "): " << refval[0] << " " << refval[1] << " " << refval[2] << " *****" << endl;
//        cout << "                 (" << ir << ", " << jr << "): " << val[0]    << " " << val[1]    << " " << val[2] << endl;
//    }
//    // END DEBUG
     
    if (!zero[0]) return false;
    if (!zero[1]) return false;
    if (!zero[2]) return false;

    return true;
}

void Contour2x2_Method_GPU::filedg4(Matrix<double> &sol_, Matrix<int> &edges_, int nedges_, 
                                int il, int jl, int ir, int jr,
                                std::vector<RealVector> &left_vrs, std::vector<RealVector> &right_vrs){

    double epsilon = 1e-10;

//    /* START_DEBUG (1/2) */
//    bool imprime = false;
//    int segmentos = 0;
//    if (nedges_ > 0) {
//        cout << "For " << nedges_ << " nedges (" << il << ", " << jl << ", " << ir << ", " << jr << ") : " << endl;
//        for(int i = 0; i < nedges_; i++){
//            cout << edges_(0, i) << " " << edges_(1, i) << " :: ";
//        }
//        cout << endl;
//        imprime = true;
//    }
//    /* END_DEBUG (1/2) The second part of DEBUG is not always necessary */

    // Store all pairs of edges that were found
    double temp[2]; temp[0] = 0.0; temp[1] = 0.0;
    RealVector p1_old(2, temp), p2_old(2, temp), p3_old(2, temp), p4_old(2, temp);
    RealVector p1(2), p2(2), p3(2), p4(2);
    for (int nedg = 0; nedg < nedges_; nedg++) {
        // LX1, LY1
        p1.component(0) = ul0 + dul * (il + sol_(0, edges_(0, nedg) ) );
        p1.component(1) = vl0 + dvl * (jl + sol_(1, edges_(0, nedg) ) );

        // LX2, LY2    
        p2.component(0) = ul0 + dul * (il + sol_(0, edges_(1, nedg) ) );
        p2.component(1) = vl0 + dvl * (jl + sol_(1, edges_(1, nedg) ) );
 
        // RX1, RY1
        p3.component(0) = ur0 + dur * (ir + sol_(2, edges_(0, nedg) ) );
        p3.component(1) = vr0 + dvr * (jr + sol_(3, edges_(0, nedg) ) );

        // RX2, RY2
        p4.component(0) = ur0 + dur * (ir + sol_(2, edges_(1, nedg) ) );
        p4.component(1) = vr0 + dvr * (jr + sol_(3, edges_(1, nedg) ) );


        /* TODO: These two "neglections" are GAMBIARRAS, HyperCube must be fixed!!! */
        // Neglect zero segments
        double norm1 = 0.0;
        double norm2 = 0.0;
        for (int comp = 0; comp < 2; comp++){
            norm1 += (p1.component(comp)-p2.component(comp))*(p1.component(comp)-p2.component(comp));
            norm2 += (p3.component(comp)-p4.component(comp))*(p3.component(comp)-p4.component(comp));
        }
        if ( (norm1 < epsilon) && (norm2 < epsilon) ) continue;

        // Neglect repetition
        norm1 = 0.0; norm2 = 0.0;
        for (int comp = 0; comp < 2; comp++){
            norm1 += (p1.component(comp)-p1_old.component(comp))*(p1.component(comp)-p1_old.component(comp))
                   + (p3.component(comp)-p3_old.component(comp))*(p3.component(comp)-p3_old.component(comp));
            norm2 += (p2.component(comp)-p2_old.component(comp))*(p2.component(comp)-p2_old.component(comp))
                   + (p4.component(comp)-p4_old.component(comp))*(p4.component(comp)-p4_old.component(comp));
        }
        if ( (norm1 < epsilon) && (norm2 < epsilon) ) continue;
        p1_old = p1; p2_old = p2; p3_old = p3; p4_old = p4;
        /* END of GAMBIARRAS!!! */

        left_vrs.push_back(p1);
        left_vrs.push_back(p2);

        right_vrs.push_back(p3);
        right_vrs.push_back(p4);

//        /* START_DEBUG (2/2) TODO: It needs the first part of the DEBUG */
//        if(imprime){
//            printf("At points (%2d, %2d, %2d, %2d) [%2d/%2d--%2d]: p1 = %1.6f, %1.6f;  p2 = %1.6f, %1.6f\n",
//                    il, jl, ir, jr, nedges_, edges_(0, nedg), edges_(1, nedg),
//                    p1.component(0), p1.component(1), p2.component(0), p2.component(1));
//            printf("                           [%2d/%2d--%2d]: p3 = %1.6f, %1.6f;  p4 = %1.6f, %1.6f\n",
//                    nedges_, edges_(0, nedg), edges_(1, nedg),
//                    p3.component(0), p3.component(1), p4.component(0), p4.component(1));
//        }
//        segmentos++;
//        /* END_DEBUG (2/2)*/

    }

//    if(nedges_ > 0) cout << "Apos gambiarras, temos " << segmentos << "/" << nedges_ << " segmentos" << endl;

    return;
}


bool Contour2x2_Method_GPU::left_right_adjacency(int il, int jl, int ir, int jr){
      return ( (fabs( (ur0+(ir+.5)*dur) - (ul0+(il+.5)*dul) ) <= dumax) &&
               (fabs( (vr0+(jr+.5)*dvr) - (vl0+(jl+.5)*dvl) ) <= dvmax) );
}


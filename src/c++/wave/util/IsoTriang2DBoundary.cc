/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) IsoTriang2DBoundary.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "IsoTriang2DBoundary.h"
using namespace std;

/*
 * ---------------------------------------------------------------
 * Definitions:
 */
using namespace std;

Three_Phase_Boundary::Three_Phase_Boundary() {
    pmin = new RealVector(2);
    pmin->component(0) = 0.0;
    pmin->component(1) = 0.0;

    pmax = new RealVector(2);
    pmax->component(0) = 1.0;
    pmax->component(1) = 1.0;

    end_edge = 1.000001;

    A_ = new RealVector(*pmin);

    B_ = new RealVector(2);

    B_->component(0) = pmin->component(0);
    B_->component(1) = pmax->component(1);


    C_ = new RealVector(2);

    C_->component(0) = pmax->component(0);
    C_->component(1) = pmin->component(1);
}

Three_Phase_Boundary::Three_Phase_Boundary(const RealVector &ppmin, const RealVector &ppmax) {
    pmin = new RealVector(ppmin);

    pmax = new RealVector(ppmax);

    end_edge = ( (pmax->component(0) + pmin->component(0)) 
               + (pmax->component(1) + pmin->component(1)) ) / 2.0 + 0.000001;

    A_ = new RealVector(*pmin);

    B_ = new RealVector(2);

    B_->component(0) = pmin->component(0);
    B_->component(1) = pmax->component(1);


    C_ = new RealVector(2);

    C_->component(0) = pmax->component(0);
    C_->component(1) = pmin->component(1);
}

Three_Phase_Boundary::Three_Phase_Boundary(const Three_Phase_Boundary &original) {
    pmin = new RealVector(*(original.pmin));

    pmax = new RealVector(*(original.pmax));

    end_edge = original.end_edge;

    A_ = new RealVector(*original.A_);

    B_ = new RealVector(*original.B_);
    C_ = new RealVector(*original.C_);
}

void Three_Phase_Boundary::envelope_curve(const FluxFunction *f, const AccumulationFunction *a,
        GridValues &gv,
        int where_constant, int number_of_steps, bool singular,
        std::vector<RealVector> &c, std::vector<RealVector> &d) {
    c.clear();
    d.clear();

    std::vector<RealVector> seg;
    edge_segments(where_constant, number_of_steps, seg);

    Envelope_Curve envelope_curve;

    envelope_curve.curve(f, a, gv, singular,
            seg,
            c, d);

    return;
}

void Three_Phase_Boundary::extension_curve(const FluxFunction *f, const AccumulationFunction *a,
        GridValues &gv,
        int where_constant, int number_of_steps, bool singular,
        int fam,int characteristic,
        std::vector<RealVector> &c, std::vector<RealVector> &d) {

    c.clear();
    d.clear();

    std::vector<RealVector> seg;
    edge_segments(where_constant, number_of_steps, seg);
    Extension_Curve extension_curve;
   
    extension_curve.curve(f, a, gv, characteristic, singular, fam,
                              seg,
                              c, d);
    
    cout<<"c: "<<c.size()<<endl;
    cout<<"d: "<<d.size()<<endl;
    return;
}

Three_Phase_Boundary::~Three_Phase_Boundary() {
    delete pmax;
    delete pmin;

    delete A_;
    delete B_;
    delete C_;
}

bool Three_Phase_Boundary::inside(const RealVector &p) const {
    return ((p.component(0) >= pmin->component(0)) && (p.component(1) >= pmin->component(1)) &&
            ((p.component(0) + p.component(1)) <= end_edge)
            );
}

bool Three_Phase_Boundary::inside(const double *p) const {
    return ((p[0] >= pmin->component(0)) && ( p[1] >= pmin->component(1)) &&
            ((p[0] + p[1]) <= end_edge)
            );
}

Boundary * Three_Phase_Boundary::clone() const {
    return new Three_Phase_Boundary(*this);
}

const RealVector& Three_Phase_Boundary::minimums(void) const {
    return *pmin;
}

const RealVector& Three_Phase_Boundary::maximums(void) const {
    return *pmax;
}

RealVector Three_Phase_Boundary::intersect(RealVector &p1, RealVector &p2) const {
    return RealVector(2);
}

const char * Three_Phase_Boundary::boundaryType() const {
    return "Three_Phase_Boundary";
}

//void Three_Phase_Boundary::physical_boundary(std::vector<RealVector> &seg){

//    std::vector<RealVector> tempSeg;

//    for (int i = 0; i < 3; i++) {
//        edge_segments(i,3,tempSeg);

//        for (int j=0; j < tempSeg.size();j++) {
//            seg.push_back(tempSeg[j]);
//        }

//    }
//}

void Three_Phase_Boundary::edge_segments(int where_constant, int number_of_steps, std::vector<RealVector> &seg) const {
    if      (where_constant == THREE_PHASE_BOUNDARY_SW_ZERO) segmented_line(*A_, *B_, number_of_steps, seg);
    else if (where_constant == THREE_PHASE_BOUNDARY_SO_ZERO) segmented_line(*A_, *C_, number_of_steps, seg);
    else                                                     segmented_line(*B_, *C_, number_of_steps, seg);

    return;
}


// THREE_PHASE_BOUNDARY_SW_ZERO 0
// THREE_PHASE_BOUNDARY_SO_ZERO 1
// THREE_PHASE_BOUNDARY_SG_ZERO 2
//
RealVector Three_Phase_Boundary::side_transverse_interior(const RealVector &p, int s) const {
    RealVector v(2);

    if (s == THREE_PHASE_BOUNDARY_SW_ZERO){
        v(0) = 1.0;
        v(1) = 0.0;
    }
    else if (s == THREE_PHASE_BOUNDARY_SO_ZERO){
        v(0) = 0.0;
        v(1) = 1.0;
    }
    else {
        v(0) = -0.707106781;
        v(1) = -0.707106781;
    }

    return v;
}

void Three_Phase_Boundary::list_of_sides(std::vector<int> &where_constant_codes, std::vector<std::string> &where_constant_names) const {
    where_constant_codes.clear();
    where_constant_names.clear();

    where_constant_codes.push_back(THREE_PHASE_BOUNDARY_SW_ZERO);
    where_constant_names.push_back(std::string("SW zero"));

    where_constant_codes.push_back(THREE_PHASE_BOUNDARY_SO_ZERO);
    where_constant_names.push_back(std::string("SO zero"));

    where_constant_codes.push_back(THREE_PHASE_BOUNDARY_SG_ZERO);
    where_constant_names.push_back(std::string("SG zero"));

    return;
}

void Three_Phase_Boundary::physical_boundary(std::vector<std::vector<RealVector> > &pb) const {
    pb.clear();

    std::vector<RealVector> side(2);

    side[0] = *A_;
    side[1] = *B_;
    pb.push_back(side);

    side[0] = *B_;
    side[1] = *C_;
    pb.push_back(side);

    side[0] = *C_;
    side[1] = *A_;
    pb.push_back(side);

    return;
}



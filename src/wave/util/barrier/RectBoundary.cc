/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RectBoundary.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include <vector>

#include "RectBoundary.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

using namespace std;

RectBoundary::RectBoundary(const RectBoundary & copy) {

    minimums_ = new RealVector(copy.minimums());
    maximums_ = new RealVector(copy.maximums());
    type_ = "rect";
    test_dimension.push_back(true);
    test_dimension.push_back(true);
    test_dimension.push_back(false); // u will not be tested for belonging to the HyperBox.

    epsilon = 0.0;
}

RectBoundary & RectBoundary::operator=(const RectBoundary & source) {

    if (this == &source)
        return *this;

    delete minimums_;
    delete maximums_;
    type_ = "rect";
    minimums_ = new RealVector(source.minimums());
    maximums_ = new RealVector(source.maximums());
    test_dimension = source.test_dimension;

    return *this;


}



int RectBoundary::edge_segments(int where_constant, int number_of_steps, const RealVector & limMin, const RealVector & limMax, std::vector<RealVector> &seg) {
    seg.clear();

    if (number_of_steps < 3) number_of_steps = 3; // Extremes must be eliminated, thus the minimum number of desired segments is 3: of which only one segment will remain.

    double p[number_of_steps + 1][3];

    double p_alpha[3], p_beta[3];

    double delta = 1.0 / (double) number_of_steps;

    //    double end_edge = 1.000001;

    // The following double is the old end_edge without the kludge.
    //    double RealEndEdge = end_edge - 0.000001;





    if (where_constant == RECT_BOUNDARY_SG_ZERO) {
        // The following protection exists for physical boundary outside the real boundary.
        //        if (pmin->component(1) > 0.0) return 0;

        p_alpha[0] = 0.0;
        p_alpha[1] = limMin(1);
        p_alpha[2] = limMin(2);
        p_beta[0] = 0.0;
        p_beta[1] = limMax(1);
        p_beta[2] = limMax(2);
    } else { // where_constant == RECT_BOUNDARY_SG_ONE
        // The following protection exists for physical boundary outside the real boundary.
        //        if (RealEndEdge < 1.0) return 0;



        p_alpha[0] = 1.0;
        p_alpha[1] = limMin(1);
        p_alpha[2] = limMin(2);
        p_beta[0] = 1.0;
        p_beta[1] = limMax(1);
        p_beta[2] = limMax(2);


    }

//
//    p_alpha[2]=1.0;
//    p_beta[2] = 1.0;


    for (int i = 0; i < number_of_steps + 1; i++) {
        double beta = (double) i*delta;
        double alpha = 1.0 - beta;
        for (int j = 0; j < 3; j++) p[i][j] = p_alpha[j] * alpha + p_beta[j] * beta;
    }

    // Boundary extension segments.
    seg.resize(2 * (number_of_steps));

    for (int i = 0; i < number_of_steps; i++) {
        seg[2 * i].resize(3);
        seg[2 * i + 1].resize(3);

        for (int j = 0; j < 3; j++) {
            seg[2 * i].component(j) = p[i][j];
            seg[2 * i + 1].component(j) = p[i + 1][j];
        }

    }

    return 1;
}

RectBoundary::RectBoundary(const RealVector & minimums, const RealVector & maximums)
: minimums_(new RealVector(minimums)),
maximums_(new RealVector(maximums)),
size_(minimums.size()), type_("rect") {
    epsilon = 0.0;
    test_dimension.push_back(true);
    test_dimension.push_back(true);
    test_dimension.push_back(false); // u will not be tested for belonging to the HyperBox.

}

RectBoundary::RectBoundary(const RealVector & minimums, const RealVector & maximums, const std::vector<bool> & test, const double eps) {
    minimums_ = new RealVector(minimums.size());
    maximums_ = new RealVector(minimums.size());

    for (int i = 0; i < minimums.size(); i++) {
        if (minimums.component(i) < maximums.component(i)) {
            minimums_->component(i) = minimums.component(i);
            maximums_->component(i) = maximums.component(i);
        } else {
            minimums_->component(i) = maximums.component(i);
            maximums_->component(i) = minimums.component(i);
        }
    }

    test_dimension.resize(minimums.size());
    for (int i = 0; i < minimums.size(); i++) test_dimension[i] = true;

    for (int i = 0; i < min((int) test.size(), minimums.size()); i++) test_dimension[i] = test[i]; //TODO ???

    //    for (int i = 0; i < test.size(); i++) test_dimension[i] = test[i];
    epsilon = eps;
}

bool RectBoundary::inside(const double *p)const {
    bool in = true;
    int pos = 0;

    while (in && pos < minimums().size()) {
        // Check if the current component should be skipped, as stated in the list of
        // exceptions.
        if (test_dimension[pos]) {
            if (p[pos] < minimums().component(pos) || p[pos] > maximums().component(pos)) in = false;
        }
        pos++;
    }

    return in;
}

bool RectBoundary::inside(const RealVector &p) const {
    bool in = true;
    int pos = 0;

    while (in && pos < minimums().size()) {
        if (test_dimension[pos]) {
            if (p.component(pos) < minimums().component(pos) || p.component(pos) > maximums().component(pos)) in = false;
        }
        //        if (p(pos) < minimums()(pos) || p(pos) > maximums()(pos)) in = false;
        pos++;
    }
    //    cout << "tamanho dentro de inside"<<in<<" "<<p.size() << endl;
    return in;


}





int RectBoundary::intersection(const RealVector &p, const RealVector &q, RealVector &r, int &w)const {


    if (inside(p) && inside(q)) return 1;
    else if (!inside(p) && !inside(q)) return -1;
    else {
        int n = minimums().size();
        double alpha, beta;
        int pos = 0;
        bool found = false;
        r.resize(p.size());

        while (pos < n && !found) {
            double d = p.component(pos) - q.component(pos);
            if (fabs(d) > epsilon * (maximums().component(pos) - minimums().component(pos))) {
                alpha = (minimums().component(pos) - q.component(pos)) / d;
                beta = (maximums().component(pos) - q.component(pos)) / d;

                if (alpha >= 0.0 && alpha <= 1.0) {
                    for (int i = 0; i < n; i++) r.component(i) = alpha * p.component(i) + (1.0 - alpha) * q.component(i);
                    found = true;
#ifdef _TEST_HYPERBOX_
                    printf("ALPHA = %f, beta = %f, pos = %d\n", alpha, beta, pos);
#endif

                    // Return the index
                    w = pos * n + 0;
                }

                if (beta >= 0.0 && beta <= 1.0) {
                    for (int i = 0; i < n; i++) r.component(i) = beta * p.component(i) + (1.0 - beta) * q.component(i);
                    found = true;
#ifdef _TEST_HYPERBOX_
                    printf("alpha = %f, BETA = %f, pos = %d\n", alpha, beta, pos);
#endif

                    // Return the index
                    w = pos * n + 1;
                }
            }
            pos++;
        }

        return 0;
    }
}

Boundary * RectBoundary::clone() const {
    return new RectBoundary(*this);
}

RectBoundary::~RectBoundary() {
    delete minimums_;
    delete maximums_;

    test_dimension.clear();

}


RealVector RectBoundary::side_transverse_interior(const RealVector &p, int s) const {
    RealVector v(2);

    if (s == RECT_BOUNDARY_SG_ZERO){
        v(0) = 1.0;
        v(1) = 0.0;
    }
    else if (s == RECT_BOUNDARY_SG_ONE){
        v(0) = -1.0;
        v(1) = 0.0;
    }

    return v;
}

void RectBoundary::edge_segments(int where_constant, int number_of_steps, std::vector<RealVector> &seg) const {
    seg.clear();

    RealVector p(2), q(2);

    if (where_constant == RECT_BOUNDARY_SG_ZERO){
        p    = *minimums_;

        q(0) = minimums_->component(0);
        q(1) = maximums_->component(1);
    }
    else if (where_constant == RECT_BOUNDARY_SG_ONE){
        p(0) = maximums_->component(0);
        p(1) = minimums_->component(1);

        q    = *maximums_;
    }
    else if (where_constant == RECT_BOUNDARY_THETA_MIN){
        p    = *minimums_;

        q(0) = maximums_->component(0);
        q(1) = minimums_->component(1);
    }
    else {
        p(0) = minimums_->component(0);
        p(1) = maximums_->component(1);

        q    = *maximums_;
    }

    segmented_line(p, q, number_of_steps, seg);

    return;
}

void RectBoundary::list_of_sides(std::vector<int> &where_constant_codes, std::vector<std::string> &where_constant_names) const {
    where_constant_codes.clear();
    where_constant_names.clear();

    where_constant_codes.push_back(RECT_BOUNDARY_SG_ZERO);
    where_constant_names.push_back(std::string("SG zero"));

    where_constant_codes.push_back(RECT_BOUNDARY_SG_ONE);
    where_constant_names.push_back(std::string("SG one"));

    where_constant_codes.push_back(RECT_BOUNDARY_THETA_MIN);
    where_constant_names.push_back(std::string("Theta min"));

    where_constant_codes.push_back(RECT_BOUNDARY_THETA_MAX);
    where_constant_names.push_back(std::string("Theta max"));

    return;
}

void RectBoundary::physical_boundary(std::vector<std::vector<RealVector> > &pb) const {
    pb.clear();

    std::vector<RealVector> temp(2);
    temp[0].resize(2);
    temp[1].resize(2);

    //
    temp[0] = *minimums_;

    temp[1](0) = minimums_->component(0);
    temp[1](1)= maximums_->component(1);

    pb.push_back(temp);

    //
    temp[0](0) = maximums_->component(0);
    temp[0](1) = minimums_->component(1);

    temp[1] = *maximums_;
    pb.push_back(temp);

    //
    temp[0] = *minimums_;

    temp[1](0) = maximums_->component(0);
    temp[1](1) = minimums_->component(1);
    pb.push_back(temp);

    //
    temp[0](0) = minimums_->component(0);
    temp[0](1) = maximums_->component(1);

    temp[1] = maximums_;
    pb.push_back(temp);

    return;
}


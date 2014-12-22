#include "Boundary.h"
#include "HugoniotContinuation.h"

Boundary::~Boundary() {
};
//double Boundary::epsilon = 1e-10;

void Boundary::segmented_line(const RealVector &p, const RealVector &q, int n, std::vector<RealVector> &seg) const{
    seg.clear();

    double delta = 1.0/(double)(n - 1);
    RealVector diff = q - p;
    
    RealVector temp = p;

    for (int i = 1; i < n; i++){
        seg.push_back(temp);

        temp = p + diff*(double)i*delta;
        seg.push_back(temp);
    }

    return;
}

int Boundary::intersection(const RealVector &p, const RealVector &q, RealVector &r, int &w) const {
    w = -1;

//    std::cout << "Boundary. p = " << p << ", q = " << q << std::endl;
//    std::cout << "inside(p) = " << inside(p) << ", inside(q) = " << inside(q) << std::endl;


    if      ( inside(p) &&  inside(q)) return BOUNDARY_INTERSECTION_BOTH_INSIDE;
    else if (!inside(p) && !inside(q)) return BOUNDARY_INTERSECTION_BOTH_OUTSIDE;
    else {
        int n = p.size();

//        // Initialize the temporary points
//        RealVector pp(p);
//        RealVector qq(q);

//        // Switch the temporary points if need be, such that pp is inside and qq is outside
//        if (!inside(pp)) {
//            RealVector temp(pp);
//            pp = qq;
//            qq = temp;
//        }

        // Point pp must be inside.
        //
        RealVector pp, qq;
        if (inside(p)){
            pp = p;
            qq = q;
        }
        else {
            pp = q;
            qq = p;
        }

        // Minimum distance.
        double d = epsilon*norm(pp - qq);

        // Iterate while the distance between the points is greater than d.
        //
        int count = 0;

        while (norm(pp - qq) > d && count < 100) {
            count++;
            r = .5*(pp + qq);

            if (inside(r)) pp = r;
            else           qq = r;
        }


        std::cout << "Boundary. Before returning, r = " << r << std::endl;
        return BOUNDARY_INTERSECTION_FOUND;

    }
}

//void Boundary::(const std::vector<RealVector> &point_on_boundary_side, int side, int n, std::vector<RealVector> &seeds){
////        if (side == THREE_PHASE_BOUNDARY_SW_ZERO){
////            // Vertical.
////            //
////            point(0) = 0.0;
////            point(1) = p;
////        }
////        else if (side == THREE_PHASE_BOUNDARY_SO_ZERO){
////            // Horizontal.
////            //
////            point(0) = p;
////            point(1) = 0.0;
////        }
////        else if (side == THREE_PHASE_BOUNDARY_SG_ZERO){
////            // Diagonal.
////            //
////            point(0) = p;
////            point(1) = 1.0 - p;
////        }

//    std::vector<double> projection;
//    for (int i = 0; i < point_on_boundary_side.size(); i++) projection.push_back(point_on_boundary_side[i](1));
//    std::sort(projection.begin(), projection.end());


//    return;
//}


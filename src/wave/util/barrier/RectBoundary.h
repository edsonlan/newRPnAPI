/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RectBoundary.h
 */

#ifndef _RectBoundary_H
#define _RectBoundary_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "Boundary.h"
#include "math.h"
#include "stdio.h"
#include <vector>
#include "Extension_Curve.h"
#include "Envelope_Curve.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

#define  RECT_BOUNDARY_SG_ZERO 0
#define  RECT_BOUNDARY_SG_ONE  1

#define  RECT_BOUNDARY_THETA_MIN 10
#define  RECT_BOUNDARY_THETA_MAX 11

class RectBoundary : public Boundary {
    private:
        RealVector * minimums_;
        RealVector * maximums_;
        int size_;
        const char * type_;

        //    std::vector<int> exception;
        std::vector<bool> test_dimension;
    public:
        RectBoundary(const RectBoundary &);
        virtual ~RectBoundary();

        Boundary * clone() const;

        RectBoundary(const RealVector & minimums, const RealVector & maximums);
        RectBoundary(const RealVector & minimums, const RealVector & maximums, const std::vector<bool> &, const double = 1e-10);

        RectBoundary & operator =(const RectBoundary &);

        double coordinateSpan(int i);
        bool inside(const RealVector &y) const;
        bool inside(const double*)const;

        const RealVector & minimums(void) const;
        const RealVector & maximums(void) const;
//        RealVector intersect(RealVector & y1, RealVector & y2) const;
    //    bool inside(const RealVector &);
    //    int intersection(const RealVector &, const RealVector &, RealVector &)const;
        int intersection(const RealVector &, const RealVector &, RealVector &, int &)const;

        const char * boundaryType()const;

        virtual int edge_segments(int where_constant, int number_of_steps, const RealVector & limMin, const RealVector & limMax, std::vector<RealVector> &seg);

       
        RealVector side_transverse_interior(const RealVector &p, int s) const;

        void edge_segments(int where_constant, int number_of_steps, std::vector<RealVector> &seg) const;
        void list_of_sides(std::vector<int> &where_constant_codes, std::vector<std::string> &where_constant_names) const;
        void physical_boundary(std::vector<std::vector<RealVector> > &pb) const;
};

inline const char * RectBoundary::boundaryType()const {
    return type_;
}



inline double RectBoundary::coordinateSpan(int i) {
    return maximums_->component(i) - minimums_->component(i);
}

//inline bool RectBoundary::inside(const RealVector & y) const {
//    // true if y inside rectangular boundary
//        cout << "tamanho dentro de inside inline" << endl;
//    bool result = true;
//    for (int i = 0; i < y.size(); i++)
//        if ((y.component(i) < minimums_->component(i)) || (y.component(i) > maximums_->component(i)))
//            result = false;
//
//    return result;
//}

inline const RealVector &RectBoundary::minimums(void) const {
    return *minimums_;
}

inline const RealVector & RectBoundary::maximums(void) const {
    return *maximums_;
}

//inline RealVector RectBoundary::intersect(RealVector & y1, RealVector & y2) const {
//    // returns a point for intersection of [y1,y2] segment with the boundary
//    double ratio;
//    RealVector vec1(size_);
//    RealVector vec2(size_);
//    RealVector result;
//
//    for (int i = 0; i < size_; i++)
//        if (y1.component(i) != y2.component(i)) {
//            ratio = (minimums_->component(i) - y1.component(i)) /      \
//                (y2.component(i) - y1.component(i));
//
//            if ((ratio >= 0) && (ratio <= 1)) {
//                //vec1.set(y1);
//                vec1 = y1;
//                //
//                // TODO
//                // RealVector doesn't have 'scale', 'add' and 'setElement' methods.
//                // In JAVA code RealVector is derivated from GVector from 'vecmath' package!
//                // daniel@impa.br
//                //
//                //vec1.scale(1d- ratio);
//                // vec2.set(y2);
//                vec2 = y2;
//                // vec2.scale(ratio);
//                // vec2.add(vec1);
//                // vec2.setElement(i, minimums_.component(i));
//                if (inside(vec2))
//                    result = vec2;
//            }
//
//            ratio = (maximums_->component(i) - y1.component(i)) /      \
//                (y2.component(i) - y1.component(i));
//
//            if ((ratio >= 0) && (ratio <= 1)) {
//                //vec1.set(y1);
//                vec1 = y1;
//                //
//                // TODO
//                // RealVector doesn't have 'scale', 'add' and 'setElement' methods.
//                // In JAVA code RealVector is derivated from GVector from 'vecmath' package!
//                // daniel@impa.br
//                //
//                //vec1.scale(1d- ratio);
//                //vec2.set(y2);
//                vec2 = y2;
//                //vec2.scale(ratio);
//                //vec2.add(vec1);
//                //vec2.setElement(i, maximums_.component(i));
//                if (inside(vec2))
//                    result = vec2;
//            }
//        }
//    return result;
//}


#endif //! _RectBoundary_H

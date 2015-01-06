/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Boundary.h
 **/

#ifndef _Boundary_H
#define	_Boundary_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */


#include "RealVector.h"
#include <iostream>
#include <vector>
#include "GridValues.h"

#define BOUNDARY_INTERSECTION_BOTH_INSIDE    1
#define BOUNDARY_INTERSECTION_FOUND          0
#define BOUNDARY_INTERSECTION_BOTH_OUTSIDE (-1)

//!

/*!
 *
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */

class Boundary {
private:
protected:
    double epsilon;

    // Divide a line whose vertices are p and q into n points.
    //
    void segmented_line(const RealVector &p, const RealVector &q, int n, std::vector<RealVector> &seg) const;
public:
    virtual ~Boundary();

    //! Check if the point is inside the boundary.
    virtual bool inside(const RealVector &y) const = 0;

    virtual bool inside(const double*) const = 0;

    //! Virtual constructor
    virtual Boundary * clone()const = 0;

    //! Minimums boundary values accessor
    virtual const RealVector & minimums() const = 0;

    //! Maximums boundary values accessor
    virtual const RealVector & maximums() const = 0;

//    virtual RealVector intersect(RealVector &y1, RealVector &y2) const = 0;

    //! Returns the boundary type
    virtual const char * boundaryType()const = 0;

    virtual int intersection(const RealVector &p, const RealVector &q, RealVector &r, int &)const;



    // Return a vector originating in p that points to the interior of the domain from 
    // the Boundary side with index s.
    //
    //    virtual RealVector side_transverse_interior(const RealVector &p, int s) const {
    //        return RealVector(2);
    //    }

    // TODO: These two methods may be unified. 
    //
    virtual RealVector side_transverse_interior(const RealVector &p, int s) const = 0;

    //    // This method must be used by all methods that need a parametrization of the sides near p.
    //    // Said parametrization must be unique.
    //    //
    //    virtual RealVector side_tangent_vector(const RealVector &p, int s) const = 0;

    double max_distance() const {
        return norm(maximums() - minimums());
    }

    // Returns the sides that form this Boundary, divided in segments.
    // To use this method the user (the GUI, for instance) needs to know 
    // what the possible values of where_constant are.
    // The method below provides those. 
    //
    virtual void edge_segments(int where_constant, int number_of_steps, std::vector<RealVector> &seg) const = 0;

    // Q.V. the commentary on the method above.
    //
    virtual void list_of_sides(std::vector<int> &where_constant_codes, std::vector<std::string> &where_constant_names) const = 0;

    // Returns the physical boundary. Derived classes may be non-regular.
    //

    virtual void physical_boundary(std::vector<std::vector<RealVector> > &pb) const {
        pb.clear();

        // Find all the sides of this Boundary.
        //
        std::vector<int> where_constant_codes;
        std::vector<std::string> where_constant_names;
        list_of_sides(where_constant_codes, where_constant_names); // Implemented in derived classes.

        for (int i = 0; i < where_constant_codes.size(); i++) {
            std::vector<RealVector> temp;
            edge_segments(where_constant_codes[i], 1000, temp); // Implemented in derived classes.

            pb.push_back(temp);
        }

        return;
    }
};

#endif	/* _Boundary_H */


/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RealSegment.h
 **/


#ifndef _RealSegment_H
#define	_RealSegment_H

#include "RealVector.h"

//!

/*! @brief Class that represents a n dimensional vector
 *
 *
 * @ingroup wave
 */


class RealSegment {
private:
    RealVector * p1_;
    RealVector * p2_;

public:

    /*! Default constructor
     */
    RealSegment(const RealVector &, const RealVector &);

    /*! Copy constructor
     */
    RealSegment(const RealSegment &);

    virtual ~RealSegment(void);


    RealSegment & operator =(const RealSegment &);
    bool operator==(const RealSegment &);


    const RealVector & p1();
    const RealVector & p2();

};

inline const RealVector & RealSegment::p1() {
    return *p1_;
}

inline const RealVector & RealSegment::p2() {
    return *p2_;
}

#endif	/* _RealSegment_H */

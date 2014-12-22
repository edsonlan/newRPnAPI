/* 
 * File:   Transform2D.h
 * Author: edsonlan
 *
 * Created on September 10, 2014, 1:52 PM
 */

#ifndef TRANSFORM2D_H
#define	TRANSFORM2D_H


#include "DoubleMatrix.h"
#include "Boundary.h"
#include <iostream>

class Transform2D {
protected:

    DoubleMatrix * transformMatrix_;
    DoubleMatrix * invertedMatrix_;

public :
    
    Transform2D(const Boundary *, int, int);
    
     Transform2D();

    virtual ~Transform2D();
    
    void virtual viewTransform( RealVector &);
    
    void virtual inverseTransform(RealVector &);

    


    void virtual viewTransform(const RealVector &, RealVector &);

    void virtual inverseTransform(const RealVector &, RealVector &);







};






#endif	/* TRANSFORM2D_H */



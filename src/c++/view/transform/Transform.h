/* 
 * File:   Transform.h
 * Author: edsonlan
 *
 * Created on November 21, 2014, 5:31 PM
 */

#ifndef TRANSFORM_H
#define	TRANSFORM_H

#include "RealVector.h"
#include "DoubleMatrix.h"


class Transform {
public:
    
    void virtual viewTransform(RealVector &)=0;

    void virtual inverseTransform(RealVector &)=0;

    
private:

};

#endif	/* TRANSFORM_H */


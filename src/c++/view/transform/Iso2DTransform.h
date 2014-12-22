/* 
 * File:   Iso2DTransform.h
 * Author: edsonlan
 *
 * Created on November 24, 2014, 11:38 AM
 */

#ifndef ISO2DTRANSFORM_H
#define	ISO2DTRANSFORM_H

#include "Transform2D.h"



class Iso2DTransform :public Transform2D{
public:
    Iso2DTransform(const Boundary &,int ,int);
    
//    
//    void viewTransform(RealVector &);
//    
//    void inverseTransform(RealVector &);
//    
//    void  viewTransform(const RealVector &, RealVector &);
//
//    void  inverseTransform(const RealVector &, RealVector &);
//
//    
//    
    
    Iso2DTransform(const Iso2DTransform& orig);
    virtual ~Iso2DTransform();
private:
    
//    
//    DoubleMatrix * transformMatrix_;
//    DoubleMatrix * invertedMatrix_;
    

};

#endif	/* ISO2DTRANSFORM_H */


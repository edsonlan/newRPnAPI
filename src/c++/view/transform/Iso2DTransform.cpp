/* 
 * File:   Iso2DTransform.cpp
 * Author: edsonlan
 * 
 * Created on November 24, 2014, 11:38 AM
 */

#include "Iso2DTransform.h"

Iso2DTransform::Iso2DTransform(const Boundary &boundary, int viewPortWidth, int viewPortHeight ):Transform2D(&boundary,viewPortWidth,viewPortHeight){
    
    
//    
//    DoubleMatrix * scaleMatrix = new DoubleMatrix(3,3);
//    
//    scaleMatrix->operator ()(0, 0) = 1.0;
//    scaleMatrix->operator ()(0, 1) = 0;
//    scaleMatrix->operator ()(0, 2) = 0;
//
//
//    scaleMatrix->operator ()(1, 0) = 0;
//    scaleMatrix->operator ()(1, 1) = sqrt(3)/2;
//    scaleMatrix->operator ()(1, 2) = 0;
//
//
//    scaleMatrix->operator ()(2, 0) = 0;
//    scaleMatrix->operator ()(2, 1) = 0;
//    scaleMatrix->operator ()(2, 2) = 1.0;
//    
//    
//    
//    DoubleMatrix * schearMatrix = new DoubleMatrix(3,3);
//    
//    schearMatrix->operator ()(0, 0) = 1.0;
//    schearMatrix->operator ()(0, 1) = 0.5;
//    schearMatrix->operator ()(0, 2) = 0;
//    
//
//
//    schearMatrix->operator ()(1, 0) = 0;
//    schearMatrix->operator ()(1, 1) = 1.0;
//    schearMatrix->operator ()(1, 2) = 0;
//
//
//    schearMatrix->operator ()(2, 0) = 0;
//    schearMatrix->operator ()(2, 1) = 0;
//    schearMatrix->operator ()(2, 2) = 1.0;
//    
//    
//    
//    
//
//
//
//DoubleMatrix tempMatrix= transformMatrix_ *(*scaleMatrix);
//
//DoubleMatrix tempMatrix2 = tempMatrix *(*schearMatrix);
//    
//    
//    delete scaleMatrix;
//    delete schearMatrix;
//    delete transformMatrix_;
//    
//    transformMatrix_ = new DoubleMatrix(tempMatrix);
    
    
//    invertedMatrix_=
    
    
    



}

Iso2DTransform::Iso2DTransform(const Iso2DTransform& orig){//:transformMatrix_(new DoubleMatrix(*orig.invertedMatrix_)),
//invertedMatrix_(new DoubleMatrix(inverse(*transformMatrix_))){
    
}


//void Iso2DTransform::viewTransform(RealVector & out){
//    
//    out = out * (*transformMatrix_);
//}
//    
//void Iso2DTransform::inverseTransform(RealVector & out){
//    
//    
//    out= out * (*invertedMatrix_);
//    
//}
//    
//void  Iso2DTransform::viewTransform(const RealVector & in, RealVector & out){
//      
//}
//
//void  Iso2DTransform::inverseTransform(const RealVector &in, RealVector & out){
//      
//      
//      
//}









Iso2DTransform::~Iso2DTransform() {
//    delete transformMatrix_;
//    delete invertedMatrix_;
}


/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RPNUMERICS.cc
 **/



//-------------------------------------
// PHYSICS
//-------------------------------------

//-------------------------------------

#include "Transform2D.h"


using namespace std;


Transform2D::Transform2D(){
    
}

Transform2D::Transform2D(const Boundary * boundary, int viewPortWidth, int viewPortHeight) : transformMatrix_(new DoubleMatrix(3, 3)) {

    //Criando transformacao

    const RealVector & boundaryMin = boundary->minimums();
    const RealVector & boundaryMax = boundary->maximums();


    double windowWidth = boundaryMax.component(0) - boundaryMin(0);

    double windowHeight = boundaryMax.component(1) - boundaryMin(1);

    double windowXOrigin = boundaryMin.component(0);
    double windowYOrigin = boundaryMin.component(1);


    double XScaleFactor = viewPortWidth / windowWidth;
    //        // not to be upside down...
    double YScaleFactor = -viewPortHeight / windowHeight;

    double XTranslateFactor = -windowXOrigin * XScaleFactor;
    // we are working with RASTER
    double YTranslateFactor = -windowYOrigin * YScaleFactor
            + viewPortHeight;
    // the viewport translation
    XTranslateFactor += 0; //?? x origem do viewPort
    YTranslateFactor += 0; // ?? y origem do viewPort



    transformMatrix_->operator()(0, 0) = XScaleFactor;
    transformMatrix_->operator()(0, 1) = 0;
    transformMatrix_->operator()(0, 2) = 0;

    transformMatrix_->operator()(1, 0) = 0;
    transformMatrix_->operator()(1, 1) = YScaleFactor;
    transformMatrix_->operator()(1, 2) = 0;

    transformMatrix_->operator()(2, 0) = XTranslateFactor;
    transformMatrix_->operator()(2, 1) = YTranslateFactor;
    transformMatrix_->operator()(2, 2) = 1;

    invertedMatrix_ = new DoubleMatrix(inverse(*transformMatrix_));



}


void Transform2D::viewTransform( RealVector & out ){
//    cout<<"transform "<< *transformMatrix_<<endl;

    out = out *(*transformMatrix_);
}

void Transform2D::inverseTransform(RealVector & out ){

//    cout<<"inverse transform "<< *invertedMatrix_<<endl;
    out = out *(*invertedMatrix_);
    
}




void Transform2D::viewTransform(const RealVector & wcCoord, RealVector & dcCoord) {


    RealVector in(3);

    in(0) = wcCoord(0);
    in(1) = wcCoord(1);
    in(2) = 1;

    RealVector outDC = in * (*transformMatrix_);
    
    
   

    dcCoord(0)=outDC(0);
    dcCoord(1)=outDC(1);

}

void Transform2D::inverseTransform(const RealVector &dcCoord, RealVector &wcCoord) {
 
    RealVector in(3);

    in(0) = dcCoord(0);
    in(1) = dcCoord(1);
    in(2) = 1;

    RealVector inWC = in*(*invertedMatrix_);
    
    wcCoord(0) = inWC(0);
    wcCoord(1) = inWC(1);

}

Transform2D::~Transform2D() {

  
    cout<<"Limpando matrizes"<<endl;
    delete transformMatrix_;
    delete invertedMatrix_;

}
/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RPNUMERICS.cc
 **/


//! Definition of RPNUMERICS.
/*!
 *
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */


#include "RpNumerics.h"





Transform2D * RpNumerics::transform_=NULL;

ViewList * RpNumerics::viewList=NULL;

int RpNumerics::curveCounter = 0;

SubPhysics * RpNumerics::physics=NULL;

 vector<RPnMethod *> * RpNumerics::methodsList_=NULL;



void RpNumerics::init(){
    
    
    physics =  new CoreyQuadSubPhysics ();
    transform_= new Iso2DTransform(*physics->boundary(),400,400);
    
    viewList= new ViewList();
    methodsList_=new vector<RPnMethod *>();
    

    
    
}


void RpNumerics::clean(){
    cout<<"Deletando"<<endl;
    
//    viewList->clean();
      delete physics;
      delete transform_;
    
      delete viewList;
      delete methodsList_;    
    
}




 Transform2D * RpNumerics::getTransform() {
    return transform_;
}
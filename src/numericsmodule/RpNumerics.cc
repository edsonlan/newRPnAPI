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


int RpNumerics::curveCounter = 0;

SubPhysics * RpNumerics::physics=NULL;

//vector <HugoniotCurve *> RpNumerics::hugoniotMethods_=NULL;



void RpNumerics::init(){
    
    
    physics =  new CoreyQuadSubPhysics ();

    
    viewList= new ViewList();
    
//    physics->list_of_Hugoniot_methods(hugoniotMethods_);
    
    
}


void RpNumerics::clean(){
    cout<<"Deletando"<<endl;
  
  
    
    
      delete physics;
      delete transform_;
    
      delete viewList;
    
    
    
    
    
    
}



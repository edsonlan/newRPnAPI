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

 ViewList * RpNumerics::viewList=NULL;
    


void RpNumerics::init(){
    
    
    physics =  new CoreyQuadSubPhysics ();
    
    viewList= new ViewList();

    
    
}


void RpNumerics::clean(){
      delete physics;
    
      delete viewList;
    
    
}



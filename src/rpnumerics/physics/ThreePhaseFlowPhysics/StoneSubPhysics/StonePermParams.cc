/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) StonePermParams.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "StonePermParams.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */




StonePermParams::StonePermParams(double expw, double expg, double expo,
                                 double expow, double expog,
                                 double cnw, double cng, double cno,
                                 double lw, double lg, 
                                 double low, double log, 
                                 double epsl):comp(new RealVector(13)) {
    
    comp->component(0)=expw;
    comp->component(1)=expg;
    comp->component(2)=expo;

    comp->component(3)=expow;
    comp->component(4)=expog;

    comp->component(5)=cnw;
    comp->component(6)=cng;
    comp->component(7)=cno;

    comp->component(8)=lw;
    comp->component(9)= lg;

    comp->component(10)= low;
    comp->component(11)= log;

    comp->component(12)= epsl;

}

StonePermParams::StonePermParams():comp(new RealVector(13)) { 
reset();
 }

StonePermParams::StonePermParams(const StonePermParams & copy):comp(new RealVector(copy.params())){
}



StonePermParams::~StonePermParams(){
    delete  comp;
}

void StonePermParams::reset() {

 comp->component(0)= 2;
    comp->component(1)=2;
    comp->component(2)=2;

    comp->component(3)=2;
    comp->component(4)=2;

    comp->component(5)=0;
    comp->component(6)=0;
    comp->component(7)=0;

    comp->component(8)=0;
    comp->component(9)= 0;

    comp->component(10)= 0;
    comp->component(11)= 0;

    comp->component(12)= 0;


}

double StonePermParams::component(int i){
    return comp->component(i);
}


/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) StoneFluxParams.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "StoneParams.h"
#include <iostream>


using namespace std;
/*
 * ---------------------------------------------------------------
 * Definitions:
 */

//const double StoneParams::DEFAULT_FLUX_PARAMS_ARRAY []= {1,1/3,1/3,1/3,0,0,0};


StoneParams::StoneParams():FluxParams(defaultParams()){}

StoneParams::StoneParams(const RealVector & params):FluxParams(params){
    grw_ = component(0);
    grg_ = component(1);
    gro_ = component(2);

    muw_ = component(3);
    mug_ = component(4);
    muo_ = component(5);

    vel_ = component(6);
}


//const FluxParams & StoneParams::DEFAULT_FLUX_PARAMS=defaultParams();






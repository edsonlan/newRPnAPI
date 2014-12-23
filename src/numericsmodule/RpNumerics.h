/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RpNumerics.h
 **/

#ifndef _RpNumerics_H
#define	_RpNumerics_H


#include "StoneSubPhysics.h" 



#include <string.h>
#include <iostream>

#include "CoreyQuadSubPhysics.h"



#include "SubPhysics.h"
#include "ViewList.h"

#include <vector>




class RpNumerics {
private:
 

    static int curveCounter;



public:
    
    static SubPhysics * physics;

    static WaveCurve * getWaveCurve(int);

    static const FluxFunction & getFlux();

    static const AccumulationFunction & getAccumulation();

    static void init();

    static void clean();

    static void clearCurveMap();

    static void removeCurve(int);
    
    static Transform2D * getTransform();

    static ViewList * viewList;
    
//    static vector<HugoniotCurve *> hugoniotMethods_;
    

};






















#endif	/* _JNIDEFS_H */


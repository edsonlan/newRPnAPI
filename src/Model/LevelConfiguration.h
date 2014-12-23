/* 
 * File:   LevelConfiguration.h
 * Author: edsonlan
 *
 * Created on December 17, 2014, 10:40 AM
 */

#ifndef LEVELCONFIGURATION_H
#define	LEVELCONFIGURATION_H

#include "Configuration.h"
#include "RealVector.h"
#include "ThreePhaseFlowEquationFunctionLevelCurve.h"
#include "RpNumerics.h"

class LevelConfiguration :public Configuration{
public:
    LevelConfiguration(const RealVector & initialPoint, int component);
    LevelConfiguration(const LevelConfiguration& orig);
    
    EquationFunctionLevelCurve * getLevelFunction()const;
   const RealVector & getInitialPoint()const ;
    int getComponent() const ;
    
    
    virtual ~LevelConfiguration();
private:

    ThreePhaseFlowEquationFunctionLevelCurve * levelCurve_;
    RealVector *  initialPoint_;
    int component_;

    
    
};





#endif	/* LEVELCONFIGURATION_H */


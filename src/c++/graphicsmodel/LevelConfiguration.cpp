/* 
 * File:   LevelConfiguration.cpp
 * Author: edsonlan
 * 
 * Created on December 17, 2014, 10:40 AM
 */

#include "LevelConfiguration.h"
#include "EquationFunctionLevelCurve.h"

LevelConfiguration::LevelConfiguration(const RealVector & initialPoint, int component): levelCurve_(new ThreePhaseFlowEquationFunctionLevelCurve(RpNumerics::physics->flux(), RpNumerics::physics->gridvalues())),
initialPoint_(new RealVector(initialPoint)), component_(component) {
}

LevelConfiguration::LevelConfiguration(const LevelConfiguration& orig) :levelCurve_(new ThreePhaseFlowEquationFunctionLevelCurve(RpNumerics::physics->flux(), RpNumerics::physics->gridvalues())),initialPoint_(new RealVector(*orig.initialPoint_)),component_(orig.component_){
    cout<<"chamando const de copia de level configuration"<<endl;
}

LevelConfiguration::~LevelConfiguration() {
    delete initialPoint_;
    delete levelCurve_;
    
}

EquationFunctionLevelCurve * LevelConfiguration::getLevelFunction()const{
    return levelCurve_;
}
const RealVector & LevelConfiguration::getInitialPoint() const{
        return *initialPoint_;
}
int LevelConfiguration::getComponent() const{
    return component_;
}
    
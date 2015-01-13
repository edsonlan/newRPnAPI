/* 
 * File:   EquationLevelSolution.cpp
 * Author: edsonlan
 * 
 * Created on December 3, 2014, 1:24 PM
 */

#include "EquationLevelMethod.h"
#include "RpNumerics.h"
#include "LevelCurveResult.h"

EquationLevelMethod::EquationLevelMethod(const LevelConfiguration & config) : RPnMethod(new LevelConfiguration(config)) {



}

RPnResult * EquationLevelMethod::calc() {

    vector<RealVector> coords;
    const LevelConfiguration & config = (const LevelConfiguration &) getConfiguration();
    config.getLevelFunction()->curve(config.getInitialPoint(), config.getComponent(), coords);
    
    Data * data = new Data(coords);
    return new LevelCurveResult(*this, data);

}

void EquationLevelMethod::recalc(Data * oldCoords) {


    
    oldCoords->clear();

    vector<RealVector> coords;
    const LevelConfiguration & config = (const LevelConfiguration &) getConfiguration();
    config.getLevelFunction()->curve(config.getInitialPoint(), config.getComponent(), coords);

    for (int j = 0; j < coords.size(); j++) {
        oldCoords->addData(coords.at(j));
        
    }


}

EquationLevelMethod::EquationLevelMethod(const EquationLevelMethod& orig) : RPnMethod(new LevelConfiguration((const LevelConfiguration &) orig.getConfiguration())) {
}

EquationLevelMethod::~EquationLevelMethod() {

}


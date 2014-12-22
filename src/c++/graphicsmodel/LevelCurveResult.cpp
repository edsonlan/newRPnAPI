/* 
 * File:   LevelCurveResult.cpp
 * Author: edsonlan
 * 
 * Created on December 18, 2014, 3:37 PM
 */

#include "LevelCurveResult.h"
#include "EquationLevelMethod.h"

LevelCurveResult::LevelCurveResult(const EquationLevelMethod & method, Data * coords):RPnResult(new EquationLevelMethod(method), coords) {
}



LevelCurveResult::~LevelCurveResult() {
}


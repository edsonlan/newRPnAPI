/* 
 * File:   LevelCurveResult.h
 * Author: edsonlan
 *
 * Created on December 18, 2014, 3:37 PM
 */

#ifndef LEVELCURVERESULT_H
#define	LEVELCURVERESULT_H

#include "RPnResult.h"
#include "EquationLevelMethod.h"

class LevelCurveResult: public RPnResult {
public:
    LevelCurveResult(const EquationLevelMethod &, Data * coords);
    virtual ~LevelCurveResult();
private:

};

#endif	/* LEVELCURVERESULT_H */


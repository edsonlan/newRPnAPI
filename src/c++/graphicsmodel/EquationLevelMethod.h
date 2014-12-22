/* 
 * File:   EquationLevelSolution.h
 * Author: edsonlan
 *
 * Created on December 3, 2014, 1:24 PM
 */

#ifndef EQUATIONLEVELSOLUTION_H
#define	EQUATIONLEVELSOLUTION_H

#include "RPnMethod.h"
#include "ThreePhaseFlowEquationFunctionLevelCurve.h"
#include "LevelConfiguration.h"

class EquationLevelMethod : public RPnMethod {
public:
    EquationLevelMethod(const LevelConfiguration &);
    EquationLevelMethod(const EquationLevelMethod& orig);


    RPnResult * calc();

   void recalc(Data *);


    virtual ~EquationLevelMethod();
private:

   

};

#endif	/* EQUATIONLEVELSOLUTION_H */


/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad4AccumulationFunction.h
 **/
#ifndef _Quad4AccumulationFunction_H
#define	_Quad4AccumulationFunction_H

//!
/*!
 *
 *
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */

#include "AccumulationFunction.h"

class Quad4AccumulationFunction : public AccumulationFunction {
    
    
public:
    
    Quad4AccumulationFunction(void);
    Quad4AccumulationFunction(const AccumulationParams & params);
    ~Quad4AccumulationFunction(void);
    
    Quad4AccumulationFunction * clone() const;
    
    int jet(const WaveState&, JetMatrix&, int) const;

};


inline Quad4AccumulationFunction::Quad4AccumulationFunction(const AccumulationParams & params) :AccumulationFunction(params) {}

#endif	/* _AccumulationFunction_H */


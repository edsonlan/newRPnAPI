/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad3AccumulationFunction.h
 **/
#ifndef _Quad3AccumulationFunction_H
#define	_Quad3AccumulationFunction_H

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

class Quad3AccumulationFunction : public AccumulationFunction {
    
    
public:
    
    Quad3AccumulationFunction(void);
    Quad3AccumulationFunction(const AccumulationParams & params);
    ~Quad3AccumulationFunction(void);
    
    Quad3AccumulationFunction * clone() const;
    
    int jet(const WaveState&, JetMatrix&, int) const;

};


inline Quad3AccumulationFunction::Quad3AccumulationFunction(const AccumulationParams & params) :AccumulationFunction(params) {}

#endif	/* _AccumulationFunction_H */


/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) AccumulationFunction.h
 **/
#ifndef _AccumulationFunction_H
#define	_AccumulationFunction_H

//!
/*!
 *
 *
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */

#include "RpFunction.h"
#include "AccumulationParams.h"

class AccumulationFunction : public RpFunction {
    
private:
    AccumulationParams * params_;
    
public:
    
    AccumulationFunction(void);
    
    AccumulationFunction(const AccumulationParams & params);
    
    AccumulationFunction(const AccumulationFunction &);
    
    virtual ~AccumulationFunction(void);
    
    void accumulationParams(const AccumulationParams & params);
    
    const AccumulationParams & accumulationParams(void) const;
    
};
inline void AccumulationFunction::accumulationParams(const AccumulationParams & params){
    
    delete params_;
    
    params_ = new AccumulationParams(params);
    
}

inline const AccumulationParams & AccumulationFunction::accumulationParams(void) const  {return *params_;}

inline AccumulationFunction::AccumulationFunction(const AccumulationParams & params) :	params_(new AccumulationParams(params)){}



#endif	/* _AccumulationFunction_H */


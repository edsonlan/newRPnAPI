/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) FluxFunction.h
 **/

#ifndef _FluxFunction_H
#define	_FluxFunction_H

#include "RpFunction.h"
#include "FluxParams.h"

//! 
/*!
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */

class FluxFunction: public RpFunction {
    
private:
    FluxParams *params_;
    
public:
   
    /*! @brief Creates a flux function with flux parameters
     * @param params  Flux parameters
     */ 

    FluxFunction();
    
    FluxFunction(const FluxParams & params);
   
    virtual ~FluxFunction(void);
    
    /*! @brief Flux parameters accessor
     *@param
     */
    
     const FluxParams & fluxParams(void) const;

    /*! @brief Flux parameters mutator 
     * @param params  New flux parameters
     */
    
    virtual void fluxParams(const FluxParams & params);
    
};




#endif	//! _FluxFunction_H

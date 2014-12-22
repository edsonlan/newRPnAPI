/* IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) FluxParams.h
 */

#ifndef _FluxParams_H
#define	_FluxParams_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RealVector.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


//!

/*!
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */


class FluxParams {
    
private:
    
    RealVector *params_;
    
    
public:
    
    /*! Constructor with size length and data array
     *@param size Length of coords array
     *@param coords Array with parameter data 
     */ 
    
    FluxParams(const int size, double *coords);

    /*! Constructor that uses a RealVector instance
     *@param params RealVector with parameter data 
     */
    
    FluxParams(const RealVector & params);
    
    /*! Copy constructor
     */
    
    FluxParams(const FluxParams & params);
    
    virtual ~FluxParams(void);
    
    /*! Parameter accessor
     */
    
    const RealVector & params(void) const ;

    /*! Parameters mutator
     *
     *@param params New parameter values 
     */
    
    void set(const RealVector & params);
    
    /*! Component parameter accessor 
     *@param index The index of parameter data array
     */
    
    double component(int index) const ;

    /*! Component parameter mutator 
     *@param index The index of parameter data array
     *@param value The new value of parameter component
     */
    
    void component(int index, double value);
    
    bool operator==(const FluxParams & fluxParams)const;
    bool operator!=(const FluxParams & fluxParams);

    FluxParams & operator=(const FluxParams &);
    
    
};

inline  const RealVector & FluxParams::params(void) const {
    return *params_;
}



inline void FluxParams::set(const RealVector & params){//TODO Create a range check
    
    delete params_;
    params_= new RealVector(params);
    
}

inline double FluxParams::component(int index) const {
    return params_->component(index);
}

inline void FluxParams::component(int index, double value) {
    params_->component(index) = value;
}

inline FluxParams & FluxParams::operator=(const FluxParams & source){
    
    if (*this== source)
        return *this;
    set(source.params());
    return *this;
    
}

inline bool FluxParams::operator!=(const FluxParams & fluxParams){
    return (!(*this==fluxParams));
}

inline bool FluxParams::operator==(const FluxParams & fluxParams) const{
    
    if (params_->size()!=fluxParams.params().size())
        return false;

    for (int i=0;i < params_->size();i++){
        
        if (params_->component(i)!=fluxParams.component(i))
            return false;
    }
    return true;
}

#endif	//! _FluxParams_H

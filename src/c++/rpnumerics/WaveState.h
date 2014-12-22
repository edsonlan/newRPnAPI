/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) WaveState.h
 */

#ifndef _WaveState_H
#define _WaveState_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */

#include "RealVector.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


/*!
 * 
 * 
 * @ingroup rpnumerics
 */

class WaveState {
private:
    
public:
    WaveState(const int dim);
    virtual ~WaveState();
    
    WaveState(const RealVector &);
    
    double operator()(const int comp) const ;
    double & operator()(const int comp);
    WaveState & operator =(const WaveState &);
    int stateSpaceDim() const ;
    void setCoords(const RealVector &);
    
protected:
    RealVector * coords_;
    
};

inline void WaveState::setCoords(const RealVector &coords){
    
    for (int i=0;i < coords_->size();i++){
        
        coords_->component(i)=coords.component(i);
    }
    
}

inline int WaveState::stateSpaceDim() const {
    
    return coords_->size();
}

inline WaveState & WaveState::operator =(const WaveState & state){
    
    if (this==&state)
        return *this;
    
    delete coords_;
    
    int i;
    
    coords_ = new RealVector(state.stateSpaceDim());
    
    for (i=0;i < coords_->size();i++){
        coords_->component(i)=state.operator ()(i);
    }
    return *this;
}




inline double WaveState::operator() (const int comp) const {
    return coords_->component(comp);
}

inline double & WaveState::operator () (const int comp) {
    return coords_->component(comp);
}



#endif //! _WaveState_H

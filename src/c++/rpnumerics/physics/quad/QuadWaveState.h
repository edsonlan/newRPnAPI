/* IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) QuadWaveState.h
 **/



//! Definition of QuadWaveState
/*!
 *
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */

#ifndef _QUAD2STATE_H
#define	_QUAD2STATE_H
#include "RealVector.h"
#include "WaveState.h"

typedef double sigma;


class QuadWaveState : public WaveState{
    
public:
    
    QuadWaveState(const int dim);
    QuadWaveState(const int dim, const sigma sig);
    QuadWaveState(const RealVector & vector);
    QuadWaveState(const RealVector & vector, const sigma sig);
        
    virtual ~QuadWaveState();

    sigma getSigma() const ;
    
    void setSigma(const sigma sig) ;
    
private:
    
    sigma sigma_;
    
    
};

inline sigma QuadWaveState::getSigma()const {return sigma_;}
inline void QuadWaveState::setSigma(const sigma sig ) { sigma_=sig;}



#endif	/* _QUAD2STATE_H */


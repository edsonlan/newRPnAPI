/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) TriphaseFluxParams.h
 */

#ifndef _StoneFluxParams_H
#define _StoneFluxParams_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */

#include "FluxParams.h"
#include "RealVector.h"
/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class StoneParams:public FluxParams {
    
private:
    double grw_, grg_, gro_;
    double muw_, mug_, muo_;
    double vel_;
    
       
    const static FluxParams & defaultParams();
public:
    
//    const static FluxParams & DEFAULT_FLUX_PARAMS;
    
    StoneParams();
    
    StoneParams(const RealVector & params);
    
    double vel()const ;
    
    double muw()const ;
    
    double muo()const ;
    
    double mug()const ;
    
    double grw()const;
    
    double gro()const ;
    
    double grg()const ;
    
};

/*
inline double StoneParams::grw()const { return params().component(0); }

inline double StoneParams::grg()const { return params().component(1); }

inline double StoneParams::gro()const { return params().component(2); }

inline double StoneParams::muw()const { return params().component(3); }

inline double StoneParams::mug()const { return params().component(4); }

inline double StoneParams::muo()const { return params().component(5); }

inline double StoneParams::vel()const { return params().component(6); }
*/

inline double StoneParams::grw()const { return grw_; }

inline double StoneParams::grg()const { return grg_; }

inline double StoneParams::gro()const { return gro_; }

inline double StoneParams::muw()const { return muw_; }

inline double StoneParams::mug()const { return mug_; }

inline double StoneParams::muo()const { return muo_; }

inline double StoneParams::vel()const { return vel_; }

inline const FluxParams & StoneParams::defaultParams()  {


    RealVector  paramsVector(7);
      
    paramsVector.component(0)=1.0;

    paramsVector.component(1)=1.0;
    paramsVector.component(2)=1.0;
    paramsVector.component(3)=1.0;
    paramsVector.component(4)=1.0;
    paramsVector.component(5)=1.0;
    paramsVector.component(6)=1.0;
    
    
    StoneParams * fluxParams = new StoneParams(paramsVector);    
    return *fluxParams;
    
}

#endif //! _StoneFluxParams_H

/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) ShockFlowParams.h
 *
 *
 **/


//! 
/*!
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */


#ifndef _ShockFlowParams_H
#define	_ShockFlowParams_H

#include "PhasePoint.h"

class ShockFlowParams {

private:
    PhasePoint * phasePoint_;
    double sigma_;
public:

    ShockFlowParams(const PhasePoint &, double);

    ShockFlowParams(const PhasePoint &);

    ShockFlowParams(const ShockFlowParams &);

    virtual ~ShockFlowParams();

    const PhasePoint & getPhasePoint() const;

    void setPhasePoint(const PhasePoint &);

    double getSigma() const;

    void setSigma(double sigma);
};





#endif	//! _ShockFlowParams_H

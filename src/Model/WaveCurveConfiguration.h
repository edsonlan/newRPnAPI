/* 
 * File:   WaveCurveConfiguration.h
 * Author: edsonlan
 *
 * Created on January 5, 2015, 4:37 PM
 */

#ifndef WAVECURVECONFIGURATION_H
#define	WAVECURVECONFIGURATION_H

#include "ReferencePoint.h"
class WaveCurveConfiguration {
public:
    WaveCurveConfiguration(const ReferencePoint &, int ,int);
    WaveCurveConfiguration(const WaveCurveConfiguration& orig);
    
    const ReferencePoint & getReferencePoint() const;
    
    int getDirection()const;
    int getFamily()const ;
    
    
    virtual ~WaveCurveConfiguration();
private:
    
    ReferencePoint * referencePoint_;
    int direction_;
    int family_;

};

#endif	/* WAVECURVECONFIGURATION_H */


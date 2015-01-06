/* 
 * File:   WaveCurveConfiguration.cpp
 * Author: edsonlan
 * 
 * Created on January 5, 2015, 4:37 PM
 */

#include "WaveCurveConfiguration.h"



WaveCurveConfiguration::WaveCurveConfiguration(const ReferencePoint & referencePoint, int family ,int direction):referencePoint_(new ReferencePoint(referencePoint)),
        family_(family),
        direction_(direction){
    
}

const ReferencePoint & WaveCurveConfiguration::getReferencePoint() const{
    return *referencePoint_;
}
    
int WaveCurveConfiguration::getDirection()const{
    return direction_;
}

int WaveCurveConfiguration::getFamily()const {
    return family_;
}



WaveCurveConfiguration::WaveCurveConfiguration(const WaveCurveConfiguration& orig):referencePoint_(new ReferencePoint(*orig.referencePoint_)),
family_(orig.family_),direction_(orig.direction_){
}

WaveCurveConfiguration::~WaveCurveConfiguration() {
    
    delete referencePoint_;
}


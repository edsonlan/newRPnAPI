/* 
 * File:   WavePointInfo.cpp
 * Author: edsonlan
 * 
 * Created on January 5, 2015, 11:27 AM
 */

#include "WavePointInfo.h"


WavePointInfo::WavePointInfo(double sigma, int eigenValueLength,double * eigenValues,Point * originPoint,MyCurve * originCurve)
:sigma_(sigma),eigenValuesLength_(eigenValueLength),eigenValues_(eigenValues),originPoint_(originPoint_),originCurve_(originCurve){
    
}





WavePointInfo::~WavePointInfo() {
    
    delete eigenValues_;
    delete originPoint_;
    delete originCurve_;
    
    
    
}


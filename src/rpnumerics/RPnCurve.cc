/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RPnCurve.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RPnCurve.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


RPnCurve::RPnCurve(const vector <RealVector> coords){
    
    resultList_=coords;
    
}
void RPnCurve::add(const RealVector  & coord){
    resultList_.push_back(coord);
}

vector<RealVector> RPnCurve::getCoords()const {
    return resultList_;
}
RPnCurve::RPnCurve(){}
int RPnCurve::size() const {
    return resultList_.size();
}

RPnCurve & RPnCurve::operator=(const RPnCurve & source) {

    resultList_.clear();
    resultList_.resize(source.getCoords().size());
    for (unsigned int i = 0; i < source.getCoords().size(); i++) {
        resultList_.at(i)=source.getCoords().at(i);
    }

    return *this;

}

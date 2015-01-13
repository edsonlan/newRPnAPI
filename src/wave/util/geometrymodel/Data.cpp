/* 
 * File:   Data.cpp
 * Author: edsonlan
 * 
 * Created on December 19, 2014, 11:04 AM
 */

#include "Data.h"

Data::Data(const vector<RealVector > & elements):elements_(new vector<RealVector>(elements)) {
    
    
    
}

Data::Data(const Data& orig) : elements_(orig.elements_) {
}

const vector<RealVector> & Data::getElements()const {
    return *elements_;
}

const RealVector & Data::getElement(int index)const {

    elements_->at(index);
}

void Data::clear() {
    elements_->clear();
}

void Data::addData(const RealVector  newElement) {
    elements_->push_back(newElement);
}

Data::~Data() {
    
//    for (std::vector<const RealVector &>::iterator it=elements_->begin(); it<elements_->end(); it++) {
//        
//        delete &*it;
//        
//        
//    }
//    
//    
    delete elements_;
    
}


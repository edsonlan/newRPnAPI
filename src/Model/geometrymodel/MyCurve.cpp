/* 
 * File:   Curve.cpp
 * Author: edsonlan
 * 
 * Created on December 19, 2014, 12:03 PM
 */

#include "MyCurve.h"

MyCurve::MyCurve(list<Data *> * data):pointsOrSegments_(data) {
}


void MyCurve::addData(Data * newElement){
pointsOrSegments_->push_back(newElement);
      
}


void MyCurve::clear(){
    
    for (std::list<Data *>::iterator it=pointsOrSegments_->begin(); it != pointsOrSegments_->end(); ++it){

      delete *it;
      
     }
    pointsOrSegments_->clear();
    
}



const Data * MyCurve::getElement(int index) const{
    
    list<Data*>::iterator it = pointsOrSegments_->begin();
    
    for (int i = 0; i < index; i++) {
        it++;

    }
    return *it;
    
        
}



void MyCurve::removeData(Data * toRemove){
    
    
    
    
}
    
const list<Data *> * MyCurve::getElements()const{
    
    return pointsOrSegments_;
}



MyCurve::MyCurve(const MyCurve& orig) {
}

MyCurve::~MyCurve() {
    
    
    cout<<"Chamando destrutor de segment"<<endl;
    
    for (std::list<Data *>::iterator it=pointsOrSegments_->begin(); it != pointsOrSegments_->end(); ++it){

      delete *it;
         
    
    }
    
    
    delete pointsOrSegments_;
    
}


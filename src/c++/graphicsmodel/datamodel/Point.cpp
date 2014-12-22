/* 
 * File:   Point.cpp
 * Author: edsonlan
 * 
 * Created on December 19, 2014, 11:14 AM
 */



#include "Point.h"

Point::Point(const RealVector & coord):coordsVector_(new list<Data *> ()),coord_(new RealVector(coord)) {
    
    coordsVector_->push_back(this);
}

const list<Data *> * Point::getElements()const{
    
    return coordsVector_;
}

const RealVector & Point::getCoord()const{
    return *coord_;
}


const Data * Point::getElement(int index) const{
        
return *coordsVector_->begin();        
        
}

void Point::clear(){
    
    delete coord_;
    
    
}


std::ostream & operator<<(std::ostream &out, const Point &r) {
    out << r.getCoord();
    return out;
}


Point::Point(const Point& orig) {
}

Point::~Point() {
    
    delete coordsVector_;
    delete coord_;
    
//    cout<<"Chamando destrutor de point"<<endl;;
    
}


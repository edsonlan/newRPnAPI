/* 
 * File:   Segment.cpp
 * Author: edsonlan
 * 
 * Created on December 19, 2014, 11:22 AM
 */

#include "Segment.h"

//Segment::Segment(const list<RealVector> & points) : points_(new list<Data *>()) {
//for (const std::list<RealVector>::iterator it=points.begin(); it != points.end(); ++it){
//
//    Point * p1 = new Point(*it);
//    points_->push_back(p1);
//
//    
//}
//
//
////    Point * p2 = new Point(points.at(1));
//
//
//
////    points_->push_back(p2);
//
//
//
//}
void Segment::clear(){
    
    for (std::list<Data *>::iterator it=points_->begin(); it != points_->end(); ++it){
        delete *it;
     }
     points_->clear();
    
}

Segment::Segment(const RealVector & p1, const RealVector & p2):points_(new list<Data *>()){
    
    
    Point * point1 = new Point(p1);
    
    Point * point2 = new Point(p2);
    
    
    points_->push_back(point1);
    points_->push_back(point2);
    
        
}

const Data * Segment::getElement(int index) const{
    
    list<Data*>::iterator it = points_->begin();
    
    for (int i = 0; i < index; i++) {
        it++;

    }
    return *it;
    
        
}

std::ostream & operator<<(std::ostream &out, const Segment &r) {
    
    
    for (std::list<Data *>::iterator it=r.points_->begin(); it != r.points_->end(); ++it){

  

     out <<*(Point *) *it;
     out<< "    ";
    
    
    }
    
    //    for (int i = 0; i < r.points_->size(); i++) {
//        Point * point =(Point *)r.points_->at(i);
//
//       
//    }
    return out;
}




const list<Data *> * Segment::getElements() const{

    return points_;
}

void Segment::removeData(Data *point) {

}

Segment::Segment(const Segment& orig) {
}

void Segment::addData(Data * point) {
    points_->push_back(point);
}

Segment::~Segment() {
//    cout<<"Chamando destrutor de segment"<<endl;
    
    for (std::list<Data *>::iterator it=points_->begin(); it != points_->end(); ++it){

      delete *it;
         
    
    }
    
    delete points_;
}
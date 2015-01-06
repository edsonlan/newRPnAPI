/* 
 * File:   Collision.cpp
 * Author: edsonlan
 * 
 * Created on January 6, 2015, 12:55 PM
 */



#include "Collision.h"
#include "CollisionException.h"

Collision::Collision():status_(0){
    
}

Collision::Collision(const RealVector & point, int status) : point_(new RealVector(point)),status_(status) {
}

Collision::Collision(const Collision& orig):point_(new RealVector(*orig.point_)),status_(orig.status_) {
}

const RealVector & Collision::getPoint()const {
    
    if (status_== 0)
        throw CollisionException();
    else
        return *point_;

}


int Collision::getStatus() const{
    return status_;
}

Collision::~Collision() {
    delete point_;
}


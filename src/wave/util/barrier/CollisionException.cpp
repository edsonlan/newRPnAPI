/* 
 * File:   CollisionException.cpp
 * Author: edsonlan
 * 
 * Created on January 6, 2015, 1:48 PM
 */

#include "CollisionException.h"

CollisionException::CollisionException() {
}

CollisionException::CollisionException(const CollisionException& orig) {
}

const char * CollisionException::what() const throw() {
    return "Collision not found ";
}
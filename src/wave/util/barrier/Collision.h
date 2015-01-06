/* 
 * File:   Collision.h
 * Author: edsonlan
 *
 * Created on January 6, 2015, 12:55 PM
 */

#ifndef COLLISION_H
#define	COLLISION_H

#include "RealVector.h"

#define INTERSECTION_FOUND          1
#define INTERSECTION_NOT_FOUND      0

class Collision {
public:
    Collision(const RealVector & ,int );
    Collision(const Collision& orig);
    Collision();
    
    int getStatus()const;
    
    const RealVector &  getPoint()const;
    
    virtual ~Collision();
private:
    
    RealVector * point_;
    
    int status_;
    
    

};

#endif	/* COLLISION_H */


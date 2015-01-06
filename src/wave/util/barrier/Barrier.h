/* 
 * File:   Barrier.h
 * Author: edsonlan
 *
 * Created on January 6, 2015, 12:55 PM
 */

#ifndef BARRIER_H
#define	BARRIER_H

#include <list>
#include "Subject.h"
#include "Segment.h"
#include "Collision.h"

using namespace std;

class Barrier:public Subject {
public:
    Barrier(list<Segment *> &);
    Barrier(const Barrier& orig);

    virtual  Collision intersect (Segment &) const =0;

    
    virtual ~Barrier();
private:

    
    list<Segment *> segments_;
    
};

#endif	/* BARRIER_H */


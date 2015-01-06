/* 
 * File:   CollisionException.h
 * Author: edsonlan
 *
 * Created on January 6, 2015, 1:48 PM
 */

#ifndef COLLISIONEXCEPTION_H
#define	COLLISIONEXCEPTION_H

#include <exception>
using namespace std;
class CollisionException:public exception {
public:
    CollisionException();
    CollisionException(const CollisionException& orig);
    
    const char * what()const throw();
//    virtual ~CollisionException();
private:

};

#endif	/* COLLISIONEXCEPTION_H */


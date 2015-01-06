/* 
 * File:   Point.h
 * Author: edsonlan
 *
 * Created on December 19, 2014, 11:14 AM
 */

#ifndef POINT_H
#define	POINT_H

#include "Data.h"


using namespace std;

class Point:public Data {
public:
    Point(const RealVector &);
    const list <Data *> * getElements() const;
    const RealVector & getCoord()const;
    
    void clear();
    
   const Data * getElement(int ) const ;
    
    friend std::ostream & operator<<(std::ostream &out, const Point &r);
    
    Point(const Point& orig);
    virtual ~Point();
private:
    
    list<Data *>  * coordsVector_;
    RealVector * coord_;;
    

};

#endif	/* POINT_H */


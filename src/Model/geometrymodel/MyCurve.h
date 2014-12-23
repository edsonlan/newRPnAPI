/* 
 * File:   Curve.h
 * Author: edsonlan
 *
 * Created on December 19, 2014, 12:03 PM
 */

#ifndef CURVE_H
#define	CURVE_H

#include "Data.h"

class MyCurve: public Data {
public:
    MyCurve(list<Data *> *);
    MyCurve(const MyCurve& orig);
    const list<Data *> * getElements()const;
    
    void clear();

    void addData(Data *);
    void removeData(Data *);
    
    
    const Data * getElement(int ) const ;
    
    
    
    virtual ~MyCurve();
private:
    
    list<Data *> * pointsOrSegments_;

};

#endif	/* CURVE_H */


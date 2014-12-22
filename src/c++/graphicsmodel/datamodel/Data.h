/* 
 * File:   Data.h
 * Author: edsonlan
 *
 * Created on December 19, 2014, 11:04 AM
 */

#ifndef DATA_H
#define	DATA_H

#include <list>
#include "RealVector.h"
using namespace std;

class Data {
public:
    Data();
    virtual const list<Data *> * getElements()const =0;
    virtual const Data * getElement(int ) const=0;
    virtual void clear()=0;
    virtual void addData(Data *){}
    virtual void removeData(Data *){}

    Data(const Data& orig);
    virtual ~Data();
private:

};

#endif	/* DATA_H */


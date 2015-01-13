/* 
 * File:   Data.h
 * Author: edsonlan
 *
 * Created on December 19, 2014, 11:04 AM
 */

#ifndef DATA_H
#define	DATA_H

#include <vector>
#include <iterator>
#include "RealVector.h"

using namespace std;

class Data {
public:
    Data(const vector<RealVector> &) ;
    const vector<RealVector> & getElements()const;
    const RealVector & getElement(int)const;
    void clear() ;
    void addData(const RealVector );


    Data(const Data& orig);
    virtual ~Data();
private:

    vector<RealVector > * elements_ ;


};

#endif	/* DATA_H */


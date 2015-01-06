/* 
 * File:   Segment.h
 * Author: edsonlan
 *
 * Created on December 19, 2014, 11:22 AM
 */

#ifndef SEGMENT_H
#define	SEGMENT_H
#include "Data.h"
#include "Point.h"

class Segment : public Data {
public:


    Segment(const RealVector &, const RealVector &);
    const Data * getElement(int)const;
    const list<Data *> * getElements()const;
    void addData(Data *);
    void removeData(Data *);
    void clear();

    friend std::ostream & operator<<(std::ostream &out, const Segment &r);


    Segment(const Segment& orig);
    virtual ~Segment();
private:

    list<Data *> * points_;


};

#endif	/* SEGMENT_H */


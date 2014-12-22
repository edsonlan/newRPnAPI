/* 
 * File:   ViewList.h
 * Author: edsonlan
 *
 * Created on September 29, 2014, 3:53 PM
 */

#ifndef VIEWLIST_H
#define	VIEWLIST_H

#include <vector>
#include "RPnMethod.h"

using namespace std;

class ViewList {
public:
    ViewList();
    ViewList(const ViewList& orig);
    void add(RPnResult *);

    void remove(RPnResult *);

    void clean();


    RPnResult * get(int);

    int getSize();

    virtual ~ViewList();
private:

    vector<RPnResult * > * curveViewList_;

};

#endif	/* VIEWLIST_H */


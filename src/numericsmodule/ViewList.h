/* 
 * File:   ViewList.h
 * Author: edsonlan
 *
 * Created on September 29, 2014, 3:53 PM
 */

#ifndef VIEWLIST_H
#define	VIEWLIST_H

#include <vector>
#include "RpnSolution.h"

using namespace std;

class ViewList {
public:
    ViewList();
    ViewList(const ViewList& orig);
    void add(RpnSolution *);

    void remove(RpnSolution *);

    void clean();


    RpnSolution * get(int);

    int getSize();

    virtual ~ViewList();
private:

    vector<RpnSolution * > * curveViewList_;

};

#endif	/* VIEWLIST_H */


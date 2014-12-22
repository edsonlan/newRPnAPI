/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RPnCurve.h
 */

#ifndef _RPnCurve_H
#define _RPnCurve_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RealVector.h"
#include <vector>
using namespace std;
/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class RPnCurve {

private:
 vector <RealVector>  resultList_;
    
public:
    
    RPnCurve ();
    RPnCurve (const vector<RealVector> );
    RPnCurve & operator=(const RPnCurve &);
    void add(const RealVector &);
    vector<RealVector> getCoords() const;
    int size() const;
    

};



#endif //! _RPnCurve_H

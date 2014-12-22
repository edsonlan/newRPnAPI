/* 
 * File:   Result.h
 * Author: edsonlan
 *
 * Created on December 16, 2014, 1:09 PM
 */

#ifndef RESULT_H
#define	RESULT_H

#include "RealVector.h"
#include "Configuration.h"
#include "Data.h"



#include <vector>
using namespace std;

class RPnMethod;

class RPnResult {
public:
    RPnResult( RPnMethod *,Data *);
    RPnResult(const RPnResult& orig);

    const Data * getData() const;
    
    virtual void realc() ;

    const RPnMethod & getMethod()const ;

    virtual ~RPnResult();
private:

    Configuration * config_;
    Data*  data_;
    RPnMethod * method_;

};

#endif	/* RESULT_H */


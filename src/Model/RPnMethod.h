/* 
 * File:   RpnSolution.h
 * Author: edsonlan
 *
 * Created on October 1, 2014, 3:31 PM
 */

#ifndef RPNMETHOD_H
#define	RPNMETHOD_H

#include "RPnResult.h"


class RPnMethod {

public:

    RPnMethod( Configuration *);
    
    virtual RPnResult * calc() {}
    
    virtual void recalc(Data *) {} //TODO so pode ser chamado pelo recalc do RpResult
    
    const Configuration & getConfiguration() const ;
    
    
    

    virtual ~RPnMethod();

private:
    
    Configuration * config_;
   

};

#endif	/* RPNSOLUTION_H */


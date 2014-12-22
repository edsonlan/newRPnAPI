/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad3FluxFunction.h
 **/

#ifndef _Quad3FluxFunction_H
#define _Quad3FluxFunction_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include <stdio.h>
#include "../../FluxFunction.h"
#include "Quad3FluxParams.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class Quad3FluxFunction : public FluxFunction {
    
public:
    Quad3FluxFunction(const Quad3FluxParams &);
    
    Quad3FluxFunction(const Quad3FluxFunction &);
    
    virtual ~Quad3FluxFunction(void);
    
    Quad3FluxFunction * clone() const ;
    
    int jet(const WaveState &u, JetMatrix &m, int degree) const;
    
    
};


#endif //! _Quad3FluxFunction_H

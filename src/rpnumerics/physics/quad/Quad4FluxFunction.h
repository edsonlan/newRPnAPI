/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad4FluxFunction.h
 **/

#ifndef _Quad4FluxFunction_H
#define _Quad4FluxFunction_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "FluxFunction.h"
#include "Quad4FluxParams.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class Quad4FluxFunction : public FluxFunction {
    
public:
    Quad4FluxFunction(const Quad4FluxParams &);
    
    Quad4FluxFunction(const Quad4FluxFunction &);
    
    virtual ~Quad4FluxFunction(void);
    
    Quad4FluxFunction * clone() const ;
    
    int jet(const WaveState &u, JetMatrix &m, int degree) const;
    
    
};


#endif //! _Quad4FluxFunction_H

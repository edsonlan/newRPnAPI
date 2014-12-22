/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad3FluxParams.h
 **/

#ifndef _Quad3FluxParams_H
#define _Quad3FluxParams_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "../../FluxParams.h"

class Quad3FluxParams : public FluxParams {

public:
	Quad3FluxParams(void);
	Quad3FluxParams(const RealVector &);
        
        virtual ~Quad3FluxParams();

//	static const double DEFAULT_A[3];
//	static const double DEFAULT_B[3][3];
//	static const double DEFAULT_C[3][3][3];

	Quad3FluxParams defaultParams(void);

};

inline Quad3FluxParams Quad3FluxParams::defaultParams(void)
{
	return Quad3FluxParams();
}


#endif //! _Quad3FluxParams_H


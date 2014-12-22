/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Quad4FluxParams.h
 **/

#ifndef _Quad4FluxParams_H
#define _Quad4FluxParams_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "FluxParams.h"

class Quad4FluxParams : public FluxParams {

public:
	Quad4FluxParams(void);
	Quad4FluxParams(const RealVector &);
        
        virtual ~Quad4FluxParams();

	static const double DEFAULT_A[2];
	static const double DEFAULT_B[2][2];
	static const double DEFAULT_C[2][2][2];

	Quad4FluxParams defaultParams(void);

};

inline Quad4FluxParams Quad4FluxParams::defaultParams(void)
{
	return Quad4FluxParams();
}


#endif //! _Quad4FluxParams_H

/*
 * JAVA CODE
 *
    public Quad2FluxParams(double[] A, double[] [] B, double[] [] [] C,int index) {
        super(Quad2.FLUX_ID, 2, A, B, C,index);
    }

    public Quad2FluxParams() {
        this(DEFAULT_A, DEFAULT_B, DEFAULT_C,0);
    }

 */

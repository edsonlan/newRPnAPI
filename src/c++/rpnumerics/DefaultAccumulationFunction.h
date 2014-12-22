/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) DefaultAccumulationFunction.h
 */

#ifndef _DefaultAccumulationFunction_H
#define _DefaultAccumulationFunction_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "AccumulationFunction.h"
#include "DefaultAccumulationParams.h"
#include "JacobianMatrix.h"
#include "HessianMatrix.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class DefaultAccumulationFunction : AccumulationFunction {

private:
	DefaultAccumulationParams accumulationParams_;

public:
	DefaultAccumulationFunction(void);

	int f(const RealVector & u, RealVector & v);
	int df(const RealVector & u, JacobianMatrix & v);
	int d2f(const RealVector & u, HessianMatrix & v);
	
	// These functions just create aliases f, df, d2f functions
	int h(const RealVector & u, RealVector & v)      { return f(u, v); }
	int dh(const RealVector & u, JacobianMatrix & v) { return df(u, v); }
	int d2h(const RealVector & u, HessianMatrix & v) { return d2f(u, v); }

	const DefaultAccumulationParams & accumulationParams(void);
	void accumulationParams(const DefaultAccumulationParams & accumulationParams);

};


inline DefaultAccumulationFunction::DefaultAccumulationFunction(void)
	: accumulationParams_(DefaultAccumulationParams())
{
}

inline int DefaultAccumulationFunction::f(const RealVector & u, RealVector & v)
{
	v = u;
	return 2;
}

inline int DefaultAccumulationFunction::df(const RealVector & u, JacobianMatrix & v)
{
	// TODO: To be implemented - daniel@impa.br
	//
	// JAVA CODE: daniel@impa.br
	//return new RealMatrix(U.size(),U.size());
	//v.equals_multiple_of_identity(1.);
	return 2;
}

inline int DefaultAccumulationFunction::d2f(const RealVector & u, HessianMatrix & v)
{
	v = * new HessianMatrix(u.size());
	return 2;
}

inline const DefaultAccumulationParams & DefaultAccumulationFunction::accumulationParams(void)
{
	return accumulationParams_;
}

inline void DefaultAccumulationFunction::accumulationParams(const DefaultAccumulationParams & accumulationParams)
{
	accumulationParams_ = accumulationParams;
}




#endif //! _DefaultAccumulationFunction_H

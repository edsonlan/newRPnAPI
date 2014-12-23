#ifndef _COREYQUAD4PHASEHUGONIOTZEROIMPLICIT_
#define _COREYQUAD4PHASEHUGONIOTZEROIMPLICIT_

#define SW_ZERO 0
#define SO_ZERO 1
#define SG_ZERO 2
#define SC_ZERO 3

#include "ZeroImplicit.h"
#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "ReferencePoint.h"

class CoreyQuad4PhaseHugoniotZeroImplicit: public ZeroImplicit {
    private:
    protected:
        RealVector (*convert)(const RealVector&);

        static RealVector sw_zero(const RealVector &p);
        static RealVector so_zero(const RealVector &p);
        static RealVector sg_zero(const RealVector &p);
        static RealVector sc_zero(const RealVector &p);

        
        
        const FluxFunction *flux;
        const AccumulationFunction *accumulation;
        ReferencePoint ref;
    public:
        CoreyQuad4PhaseHugoniotZeroImplicit();
        virtual ~CoreyQuad4PhaseHugoniotZeroImplicit();

        RealVector f_zero(const RealVector &p);
        void find_zeroes(int side, std::vector<RealVector> &zeroes);
};

#endif // _COREYQUAD4PHASEHUGONIOTZEROIMPLICIT_


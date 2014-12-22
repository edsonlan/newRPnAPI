#ifndef _SECONDARY_BIFURCATION_
#define _SECONDARY_BIFURCATION_

#include "ThreeImplicitFunctions.h"
#include "Contour2x2_Method.h"

class Secondary_Bifurcation : public ThreeImplicitFunctions {
    private:
    protected:
        Matrix<double> flux_left, accum_left;

        const FluxFunction *lff;
    	const AccumulationFunction *laa;

        const FluxFunction *rff;
        const AccumulationFunction *raa;

    public:
        Secondary_Bifurcation(){
            gv_left = gv_right = 0;

            flux_left.resize(2, 4);
            accum_left.resize(2, 4);

            singular = true;
        }

        ~Secondary_Bifurcation(){}

        bool prepare_cell(int i, int j);

        bool function_on_cell(double *val, int ir, int jr, int kl, int kr);

        void curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues &lg,
                   const FluxFunction *rf, const AccumulationFunction *ra, GridValues &rg,
                   std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve);
};

#endif // _SECONDARY_BIFURCATION_

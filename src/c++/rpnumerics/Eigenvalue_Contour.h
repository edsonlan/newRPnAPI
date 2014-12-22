#ifndef _EIGENVALUE_CONTOUR_
#define _EIGENVALUE_CONTOUR_

#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Matrix.h"
#include "RealVector.h"
#include "eigen.h"

#include "ImplicitFunction.h"
#include "ContourMethod.h"

#include <algorithm>

class Eigenvalue_Contour : public ImplicitFunction {
    private:
    protected:
        std::vector<double> levels;
        double level;
        int family;

        // Find the minimum and maximum lambdas, as were computed on the grid.
        //
        void find_minmax_lambdas(int family, double &min, double &max);
    public:
        Eigenvalue_Contour();
        ~Eigenvalue_Contour(){}

        // Set the levels as a vector with arbitrary values.
        // The levels will be sorted from minimum to maximum.
        //
        void set_levels(int family, const std::vector<double> &l);

        // Set the levels from the minimum towards the maximum, with a
        // uniform separation between levels.
        //
        void set_levels_from_delta(int family, double delta_l);

        // Set the number of levels to be uniformly distributed.
        //
        void set_number_levels(int family, int n);

        // Set the level for a given family.
        //
        void set_level(double l, int family);

        // Set the level at the given point.
        //
        void set_level_from_point(const FluxFunction *f, const AccumulationFunction *a,
                                  int family, const RealVector &p);

        // Set the levels radiating from the level at the given point, 
        // with the given distance between levels.
        //
        void set_levels_from_point(const FluxFunction *f, const AccumulationFunction *a,
                                   GridValues &g, 
                                   int family, const RealVector &p, double delta_l);

        int function_on_square(double *foncub, int i, int j);

        int curve(const FluxFunction *f, const AccumulationFunction *a, 
                  GridValues &g, 
                  std::vector< std::vector<RealVector> > &curve, std::vector<double> &out_levels);
        int curve(const FluxFunction *f, const AccumulationFunction *a, 
                  GridValues &g,
                  std::vector<RealVector> &curve, double &l);
};

#endif // _EIGENVALUE_CONTOUR_


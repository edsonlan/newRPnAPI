#ifndef _COINCIDENCE_CONTOUR_
#define _COINCIDENCE_CONTOUR_

#include "Coincidence.h"
#include "ImplicitFunction.h"
#include "ContourMethod.h"

#define CHARACTERISTIC_SPEED_EVAPORATION 0
#define CHARACTERISTIC_SPEED_SATURATION  1

class Coincidence_Contour : public ImplicitFunction {
    private:
    protected:
        const Coincidence *coincidence;

        static int coincidence_on_square(Coincidence_Contour *obj, double *foncub, int i, int j);

        static int evaporation_on_square(Coincidence_Contour *obj, double *foncub, int i, int j);
        double evaporation_level; // To be improved as a grid.

        static int saturation_on_square(Coincidence_Contour *obj, double *foncub, int i, int j);
        double saturation_level; // To be improved as a grid.

        int (*fos)(Coincidence_Contour *obj, double *foncub, int i, int j);
    public:
        Coincidence_Contour(const Coincidence *c);
        ~Coincidence_Contour();

        int function_on_square(double *foncub, int i, int j){
            return (*fos)(this, foncub, i, j);
        }

        int curve(const FluxFunction *f, const AccumulationFunction *a, 
                  GridValues &g, std::vector<RealVector> &coincidence_curve);

        int characteristic_speed_curve(const FluxFunction *f, const AccumulationFunction *a, 
                                       GridValues &g, 
                                       const RealVector &p, int type, 
                                       std::vector<RealVector> &curve, double &lev);
};

#endif // _COINCIDENCE_CONTOUR_


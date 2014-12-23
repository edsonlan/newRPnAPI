#ifndef _QUAD2EXPLICITHUGONIOTCURVE_
#define _QUAD2EXPLICITHUGONIOTCURVE_

#include "HugoniotCurve.h"
#include "ParametricPlot.h"
#include "Utilities.h"

#define QUAD2_GENERIC_POINT 0

class Quad2SubPhysics;

class Quad2ExplicitHugoniotCurve : public HugoniotCurve {
    private:
    protected:
        static double sign(double v){
            return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0.0);
        }

        static RealVector generic(void *cqeh, double phi);
        Quad2SubPhysics *quad2;

        double alpha0, alpha1, alpha2;
        double beta0,  beta1,  beta2;
        double gamma0, gamma1, gamma2;
    public:
        Quad2ExplicitHugoniotCurve(Quad2SubPhysics *c);
        virtual ~Quad2ExplicitHugoniotCurve();

        void curve(const ReferencePoint &ref, int type, std::vector<Curve> &c);

        static bool f_asymptote(void *obj, const RealVector &p, const RealVector &q);

        void list_of_reference_points(std::vector<int> &type, std::vector<std::string> &name) const {
            type.clear();
            type.push_back(QUAD2_GENERIC_POINT);

            name.clear();
            name.push_back(std::string("Generic point"));

            return;
        }
};

#endif // _QUAD2EXPLICITHUGONIOTCURVE_


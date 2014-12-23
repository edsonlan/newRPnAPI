#ifndef _EXPLICITTRANSITIONALWAVECURVE_
#define _EXPLICITTRANSITIONALWAVECURVE_

#include "CoreyQuadSubPhysics.h"
#include "Bisection.h"
#include "Utilities.h"
#include "eigen.h"

#define EW 0
#define BO 1
#define DG 2

#define CURVE_OK      0
#define CURVE_ERROR (-1)

#define SUBDIVISION_OK      0
#define SUBDIVISION_ERROR (-1)

class ExplicitTransitionalWavecurve {
    private:
    protected:
        CoreyQuadSubPhysics *corey_;

        // This member is needed by the method below.
        //
        double r, sM;

        void transitional_double_contact_jet(double L, int degree, JetMatrix &Ljet);
        void U_side_polynomial_jet(double L, int degree, JetMatrix &Ljet);

        RealVector EW_line(double sM);
        double EW_s(const RealVector &M);

        RealVector BO_line(double sM);
        double BO_s(const RealVector &M);

        RealVector DG_line(double sM);
        double DG_s(const RealVector &M);
    public:
        ExplicitTransitionalWavecurve();
        virtual ~ExplicitTransitionalWavecurve();

        static int L_function(void *obj, double L, double &C2);
        static int U_side_function(void *obj, double L, double &sl);

        int subdivide_curve(int type, std::vector<double> &s_param, std::vector<RealVector> &extpts);
        int find_subdivision(int type, const RealVector &M, const std::vector<double> &s_param, RealVector &point_sl, RealVector &point_sr);
};

#endif // _EXPLICITTRANSITIONALWAVECURVE_


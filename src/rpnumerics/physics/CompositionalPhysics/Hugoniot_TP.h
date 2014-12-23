#ifndef _HUGONIOT_TP_
#define _HUGONIOT_TP_

#include "ReferencePoint.h"
//#include "Hugoniot_Locus.h"

#include "SubPhysics.h"
#include "HugoniotCurve.h"
#include "ImplicitFunction.h"

#define HUGONIOTTP_GENERIC_POINT 0

class Hugoniot_TP : public HugoniotCurve, public ImplicitFunction {
    private:
    protected:
        RealVector Fref, Gref;
        RealVector Uref;

        DoubleMatrix JFref, JGref;

        const FluxFunction *ff;
        const AccumulationFunction *aa;

        SubPhysics *subphysics_;
    public:
        Hugoniot_TP(const FluxFunction *ff, const AccumulationFunction *aa, const Boundary *bb);
        ~Hugoniot_TP();

        void subphysics(SubPhysics *s){
            subphysics_ = s;

            return;
        }

        double complete_points(const RealVector &Uplus);

        int function_on_square(double *foncub, int i, int j);

        // TODO: To be deprecated.
        //
        int classified_curve(const FluxFunction *f, const AccumulationFunction *a, 
                             GridValues &g, /* const RealVector &r, */
                             const ReferencePoint &r,
                             std::vector<HugoniotPolyLine> &hugoniot_curve);

        // TODO: To be deprecated.
        //
        int classified_curve(const FluxFunction *f, const AccumulationFunction *a, 
                             GridValues &g, /* const RealVector &r, */
                             const ReferencePoint &r,
                             std::vector<HugoniotPolyLine> &hugoniot_curve,
                             std::vector<bool> &circular);
        

        // TODO: To be deprecated.
        //
        int curve(const FluxFunction *f, const AccumulationFunction *a, 
                  GridValues &g, 
                  const RealVector &r,
                  std::vector<RealVector> &hugoniot_curve);

        // TODO: To be deprecated.
        //
        int curve(const FluxFunction *f, const AccumulationFunction *a, 
                  GridValues &g, 
/*                  const FluxFunction *ref_f, const AccumulationFunction *ref_a,
                  const RealVector &r,*/
                  const ReferencePoint &r,
                  std::vector<RealVector> &hugoniot_curve);

        void curve(const ReferencePoint &ref, int type, std::vector<Curve> &c);

        void map(const RealVector &p, double &f, RealVector &map_Jacobian);

        bool improvable(void);

        virtual void list_of_reference_points(std::vector<int> &type, std::vector<std::string> &name) const {
            type.clear();
            name.clear();

            type.push_back(HUGONIOTTP_GENERIC_POINT);
            name.push_back(std::string("Generic point"));

            return;
        }
};

#endif // _HUGONIOT_TP_


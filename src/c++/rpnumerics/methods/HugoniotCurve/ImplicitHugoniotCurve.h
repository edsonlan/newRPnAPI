#ifndef _IMPLICITHUGONIOTCURVE_
#define _IMPLICITHUGONIOTCURVE_

#include "HugoniotCurve.h"
#include "ImplicitFunction.h"
#include "ContourMethod.h"
#include "Newton_Improvement.h"

// Implicit Hugoniots need a SubPhysics that will provide them
// with a working GridValues.
//
#include "SubPhysics.h"

#define IMPLICITHUGONIOTCURVE_GENERIC_POINT 0

class ImplicitHugoniotCurve : public HugoniotCurve, public ImplicitFunction {
    private:
    protected:
        // TODO: Maybe these four members below can be eliminated.
        //
        RealVector Fref, Gref;
        DoubleMatrix JFref, JGref;

        SubPhysics *subphysics_;
    public:
        ImplicitHugoniotCurve(const FluxFunction *ff, const AccumulationFunction *aa, const Boundary *bb);
        virtual ~ImplicitHugoniotCurve();

        void subphysics(SubPhysics *s){
            subphysics_ = s;

            return;
        }

        SubPhysics* subphysics(){
            return subphysics_;
        }

        // TODO: To be deprecated.
        //
        void set_grid(GridValues *g){gv = g; return;}
    
        void curve(const ReferencePoint &ref, int type, std::vector<Curve> &c);

        virtual int function_on_square(double *foncub, int i, int j);
        bool improvable(void);
        int complete(const RealVector &p0, const RealVector &p1, const RealVector &p_init, RealVector &p_completed);

        // TODO: Eliminate this method.
        //
        void types(std::vector<int> &type, std::vector<std::string> &name) const {
            type.clear();
            name.clear();

            return;
        }

        virtual void list_of_reference_points(std::vector<int> &type, std::vector<std::string> &name) const {
            type.clear();
            name.clear();

            type.push_back(IMPLICITHUGONIOTCURVE_GENERIC_POINT);
            name.push_back(std::string("Generic point"));

            return;
        }

        void curve(const ReferencePoint &ref, int type, std::vector<HugoniotPolyLine> &classified_curve){
            HugoniotCurve::curve(ref, type, classified_curve);

            return;
        }

        double speed(const RealVector &F0, const RealVector &G0, const RealVector &F1, const RealVector &G1);
        void equilibria(const std::vector<Curve> &curve, const ReferencePoint& ref, const RealVector &p, std::vector<RealVector> &ep);

        // See if this can be improved.
        //
        GridValues* gridvalues(){return gv;}

        RealVector& F_ref(){return Fref;}
        RealVector& G_ref(){return Gref;}

        DoubleMatrix& JF_ref(){return JFref;}
        DoubleMatrix& JG_ref(){return JGref;}
};

#endif // _IMPLICITHUGONIOTCURVE_


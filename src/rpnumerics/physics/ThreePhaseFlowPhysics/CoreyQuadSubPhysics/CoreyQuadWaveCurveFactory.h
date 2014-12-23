#ifndef _COREYQUADWAVECURVEFACTORY_
#define _COREYQUADWAVECURVEFACTORY_

#include "WaveCurveFactory.h"
#include "CoreyQuadSubPhysics.h"

#define COREYQUADWAVECURVEFACTORY_GENERIC_POINT WAVECURVEFACTORY_GENERIC_POINT

#define COREYQUADWAVECURVEFACTORY_O_TO_B        1001
#define COREYQUADWAVECURVEFACTORY_O_TO_W        1002
#define COREYQUADWAVECURVEFACTORY_O_TO_G        1003

#define COREYQUADWAVECURVEFACTORY_W_TO_E        1004
#define COREYQUADWAVECURVEFACTORY_W_TO_G        1005
#define COREYQUADWAVECURVEFACTORY_W_TO_O        1006

#define COREYQUADWAVECURVEFACTORY_G_TO_D        1007
#define COREYQUADWAVECURVEFACTORY_G_TO_W        1008
#define COREYQUADWAVECURVEFACTORY_G_TO_O        1009

#define COREYQUADWAVECURVEFACTORY_GW_SIDE       1010
#define COREYQUADWAVECURVEFACTORY_GO_SIDE       1011
#define COREYQUADWAVECURVEFACTORY_WO_SIDE       1012

#define COREYQUADWAVECURVEFACTORY_INVALID_PARAMETERS WAVECURVEFACTORY_INVALID_PARAMETERS

class CoreyQuadWaveCurveFactory: public WaveCurveFactory {
    private:
    protected:
        CoreyQuadSubPhysics *coreyquadsubphysics_;
    public:
        CoreyQuadWaveCurveFactory(const FluxFunction *ff, const AccumulationFunction *gg, 
                                       const Boundary *bb, const ODE_Solver *o, 
                                       RarefactionCurve *r, ShockCurve *s, CompositeCurve *c, 
                                       CoreyQuadSubPhysics *tpfsp);
        virtual ~CoreyQuadWaveCurveFactory();

        virtual int wavecurve(int type, const RealVector &initial_point, int family, int increase, HugoniotContinuation *h, WaveCurve &hwc, 
                              int &wavecurve_stopped_because, int &edge){

            return wavecurve(type, initial_point, family, increase, h, 
                             0, 0, 
                             hwc, 
                             wavecurve_stopped_because, edge);

        }

        virtual int wavecurve(int type, const RealVector &initial_point, int family, int increase, HugoniotContinuation *h, 
                              void *linobj, double (*linear_function)(void *o, const RealVector &p),
                              WaveCurve &hwc, 
                              int &wavecurve_stopped_because, int &edge);

        void list_of_initial_points(std::vector<int> &type, std::vector<std::string> &name) const;
};

#endif // _COREYQUADWAVECURVEFACTORY_


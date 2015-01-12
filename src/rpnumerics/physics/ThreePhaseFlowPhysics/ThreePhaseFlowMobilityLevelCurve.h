#ifndef _THREEPHASEFLOWMOBILITYLEVELCURVE_
#define _THREEPHASEFLOWMOBILITYLEVELCURVE_

#include "ImplicitFunction.h"
#include "ThreePhaseFlowMobility.h"
#include "ThreePhaseFlowSubPhysics.h"

#include "ContourMethodPure.h"

#define WATER_MOBILITY_CURVE 0
#define OIL_MOBILITY_CURVE   1
#define GAS_MOBILITY_CURVE   2

class ThreePhaseFlowMobilityLevelCurve: public ImplicitFunction {
    private:
    protected:
        ThreePhaseFlowMobility   *mobility_;
        ThreePhaseFlowSubPhysics *subphysics_;

        static double water_mobility(ThreePhaseFlowMobilityLevelCurve *obj, const RealVector &p);
        static double oil_mobility(ThreePhaseFlowMobilityLevelCurve *obj, const RealVector &p);
        static double gas_mobility(ThreePhaseFlowMobilityLevelCurve *obj, const RealVector &p);

        double (*mobilityfunction)(ThreePhaseFlowMobilityLevelCurve *obj, const RealVector &p);
        double level_;
        int component_;

        Matrix<RealVector> mobility_on_grid;
        void init();
    public:
        ThreePhaseFlowMobilityLevelCurve(ThreePhaseFlowSubPhysics *s);
        virtual ~ThreePhaseFlowMobilityLevelCurve();

        virtual int function_on_square(double *foncub, int i, int j);
        void curve(const RealVector &ref, int type, std::vector<RealVector> &c);
        void curve(double level, int type, std::vector<RealVector> &c);

        double level(const RealVector &ref, int type);

        void list_of_types(std::vector<int> &type, std::vector<std::string> &name){
            type.clear();
            type.push_back(WATER_MOBILITY_CURVE);
            type.push_back(OIL_MOBILITY_CURVE);
            type.push_back(GAS_MOBILITY_CURVE);

            name.clear();
            name.push_back(std::string("Water mobility"));
            name.push_back(std::string("Oil mobility"));
            name.push_back(std::string("Gas mobility"));
        }

        void prepare_grid(){init(); return;}
};

#endif // _THREEPHASEFLOWMOBILITYLEVELCURVE_


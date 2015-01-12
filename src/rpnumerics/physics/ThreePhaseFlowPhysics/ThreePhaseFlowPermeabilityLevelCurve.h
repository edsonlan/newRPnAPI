#ifndef _THREEPHASEFLOWPERMEABILITYLEVELCURVE_
#define _THREEPHASEFLOWPERMEABILITYLEVELCURVE_

#include "ImplicitFunction.h"
#include "ThreePhaseFlowPermeability.h"
#include "ThreePhaseFlowSubPhysics.h"

#include "ContourMethodPure.h"

#define WATER_PERMEABILITY_CURVE 0
#define OIL_PERMEABILITY_CURVE   1
#define GAS_PERMEABILITY_CURVE   2

class ThreePhaseFlowPermeabilityLevelCurve: public ImplicitFunction {
    private:
    protected:
        ThreePhaseFlowPermeability *permeability_;
        ThreePhaseFlowSubPhysics   *subphysics_;

        static double water_permeability(ThreePhaseFlowPermeabilityLevelCurve *obj, const RealVector &p);
        static double oil_permeability(ThreePhaseFlowPermeabilityLevelCurve *obj, const RealVector &p);
        static double gas_permeability(ThreePhaseFlowPermeabilityLevelCurve *obj, const RealVector &p);

        double (*permeabilityfunction)(ThreePhaseFlowPermeabilityLevelCurve *obj, const RealVector &p);
        double level_;
        int component_;

        Matrix<RealVector> permeability;
        void init();
    public:
        ThreePhaseFlowPermeabilityLevelCurve(ThreePhaseFlowSubPhysics *s);
        ThreePhaseFlowPermeabilityLevelCurve(ThreePhaseFlowSubPhysics *s, ThreePhaseFlowPermeability *p);
        virtual ~ThreePhaseFlowPermeabilityLevelCurve();

        virtual int function_on_square(double *foncub, int i, int j);
        void curve(const RealVector &ref, int type, std::vector<RealVector> &c);
        void curve(double level, int type, std::vector<RealVector> &c);

        double level(const RealVector &ref, int type);

        void list_of_types(std::vector<int> &type, std::vector<std::string> &name){
            type.clear();
            type.push_back(WATER_PERMEABILITY_CURVE);
            type.push_back(OIL_PERMEABILITY_CURVE);
            type.push_back(GAS_PERMEABILITY_CURVE);

            name.clear();
            name.push_back(std::string("Water permeability"));
            name.push_back(std::string("Oil permeability"));
            name.push_back(std::string("Gas permeability"));
        }

        void prepare_grid(){init(); return;}
};

#endif // _THREEPHASEFLOWPERMEABILITYLEVELCURVE_


#ifndef _THREEPHASEFLOWEQUATIONFUNCTIONLEVELCURVE_
#define _THREEPHASEFLOWEQUATIONFUNCTIONLEVELCURVE_

#include "EquationFunctionLevelCurve.h"

#define WATER_COMPONENT 0
#define OIL_COMPONENT 1
#define GAS_COMPONENT 2

class ThreePhaseFlowEquationFunctionLevelCurve: public EquationFunctionLevelCurve {
    private:
    protected:
        static double fw_level_function(EquationFunctionLevelCurve *obj, const RealVector &p);
        static double fo_level_function(EquationFunctionLevelCurve *obj, const RealVector &p);
        static double fg_level_function(EquationFunctionLevelCurve *obj, const RealVector &p);
    public:
        ThreePhaseFlowEquationFunctionLevelCurve(const RpFunction *rpf, GridValues *g);
        virtual ~ThreePhaseFlowEquationFunctionLevelCurve();

        virtual void curve(const RealVector &ref, int component, std::vector<RealVector> &c);

        virtual double level(const RealVector &ref, int component);

        void list_of_components(std::vector<int> &component, std::vector<std::string> &name){
            component.clear();
            component.push_back(WATER_COMPONENT);
            component.push_back(OIL_COMPONENT);
            component.push_back(GAS_COMPONENT);

            name.clear();
            name.push_back(std::string("Water component"));
            name.push_back(std::string("Oil component"));
            name.push_back(std::string("Gas component"));

        }
};

#endif // _THREEPHASEFLOWEQUATIONFUNCTIONLEVELCURVE_


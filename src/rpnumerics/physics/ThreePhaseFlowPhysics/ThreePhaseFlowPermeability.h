#ifndef _THREEPHASEFLOWPERMEABILITY_
#define _THREEPHASEFLOWPERMEABILITY_

#include "AuxiliaryFunction.h"

#define JETTESTER_ENABLED_PERMEABILITY


class ThreePhaseFlowSubPhysics;

class ThreePhaseFlowPermeability: public AuxiliaryFunction 
                                
{
    private:
    protected:
        ThreePhaseFlowSubPhysics *subphysics_;
    public:
        ThreePhaseFlowPermeability(ThreePhaseFlowSubPhysics *s);
        virtual ~ThreePhaseFlowPermeability();

        virtual int PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water) = 0;
        virtual int PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas) = 0;
        virtual int PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil) = 0;

        virtual void reduced_permeability(const RealVector &state, RealVector &rp) = 0;

        #ifdef JETTESTER_ENABLED_PERMEABILITY
        static int testable_PermeabilityWater_jet(void *obj, const RealVector &state, int degree, JetMatrix &jm){
            ThreePhaseFlowPermeability *tpfp = (ThreePhaseFlowPermeability*)obj;
            return tpfp->PermeabilityWater_jet(state, degree, jm);
        }

        static int testable_PermeabilityOil_jet(void *obj, const RealVector &state, int degree, JetMatrix &jm){
            ThreePhaseFlowPermeability *tpfp = (ThreePhaseFlowPermeability*)obj;
            return tpfp->PermeabilityOil_jet(state, degree, jm);
        }

        static int testable_PermeabilityGas_jet(void *obj, const RealVector &state, int degree, JetMatrix &jm){
            ThreePhaseFlowPermeability *tpfp = (ThreePhaseFlowPermeability*)obj;
            return tpfp->PermeabilityGas_jet(state, degree, jm);
        }

        void list_of_functions(std::vector<int (*)(void*, const RealVector&, int degree, JetMatrix&)> &list,
                               std::vector<std::string> &name){
            list.clear();

            list.push_back(&testable_PermeabilityWater_jet);
            list.push_back(&testable_PermeabilityOil_jet);
            list.push_back(&testable_PermeabilityGas_jet);

            name.clear();
            name.push_back(std::string("Water permeability"));
            name.push_back(std::string("Oil permeability"));
            name.push_back(std::string("Gas permeability"));

            return;
        }
        #endif  
};

#endif // _THREEPHASEFLOWPERMEABILITY_


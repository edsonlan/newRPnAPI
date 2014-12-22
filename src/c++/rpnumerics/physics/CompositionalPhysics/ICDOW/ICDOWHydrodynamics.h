#ifndef _ICDOWHYDRODYNAMICS_
#define _ICDOWHYDRODYNAMICS_

#include "Parameter.h"
#include "JetMatrix.h"

#define JETTESTER_ENABLED_ICDOWHYDRO

#ifdef JETTESTER_ENABLED_ICDOWHYDRO
//#include "JetTester.h"

#endif

class ICDOWHydrodynamics 
                        
{
    private:
    protected:
        Parameter *muw_parameter, *muo_parameter;
        Parameter *grw_parameter, *gro_parameter;
        Parameter *vel_parameter;
        Parameter *swc_parameter, *lambda_parameter;
    public:
        ICDOWHydrodynamics(Parameter *muw, Parameter *muo,
                           Parameter *grw, Parameter *gro,
                           Parameter *vel,
                           Parameter *swc, Parameter *lambda);
        virtual ~ICDOWHydrodynamics();

        void water_fractional_flow(double sw, int degree, JetMatrix &fw_jet);

        void oil_fractional_flow(double sw, int degree, JetMatrix &fo_jet);

        #ifdef JETTESTER_ENABLED_ICDOWHYDRO
        static int testable_water_jet(void *obj, const RealVector &state, int degree, JetMatrix &jm){
            ICDOWHydrodynamics *icdowhydro = (ICDOWHydrodynamics*)obj;

            double sw = state(0);
            icdowhydro->water_fractional_flow(sw, degree, jm);

            return 0;
        }

        static int testable_oil_jet(void *obj, const RealVector &state, int degree, JetMatrix &jm){
            ICDOWHydrodynamics *icdowhydro = (ICDOWHydrodynamics*)obj;

            double sw = state(0);
            icdowhydro->oil_fractional_flow(sw, degree, jm);

            return 0;
        }

        void list_of_functions(std::vector<int (*)(void*, const RealVector&, int degree, JetMatrix&)> &list,
                               std::vector<std::string> &name){
            list.clear();
            list.push_back(&testable_water_jet);
            list.push_back(&testable_oil_jet);

            name.clear();
            name.push_back(std::string("ICDOW Hydrodynamics (water)"));
            name.push_back(std::string("ICDOW Hydrodynamics (oil)"));

            return;
        }
        #endif
};

#endif // _ICDOWHYDRODYNAMICS_


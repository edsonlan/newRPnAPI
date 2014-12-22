#ifndef _DEADVOLATILEVOLATILEGASTHERMODYNAMICS_
#define _DEADVOLATILEVOLATILEGASTHERMODYNAMICS_

#include "AuxiliaryFunction.h"

class DeadVolatileVolatileGasThermodynamics : public AuxiliaryFunction {
    private:
    protected:
        Parameter *B_parameter, *D_parameter;
        Parameter *mu_oB_parameter, *mu_oD_parameter, *mu_G_parameter;
        Parameter *rg_parameter, *re_parameter;
    public:
        DeadVolatileVolatileGasThermodynamics(Parameter *B, Parameter *D, Parameter *mu_oB, Parameter *mu_oD, Parameter *mu_G, Parameter *rg, Parameter *re);
        virtual ~DeadVolatileVolatileGasThermodynamics();

        void viscosity_ratio(int degree, double y, JetMatrix &R);

        void gas_molar_density_a(int degree, double y, JetMatrix &Rga_jet);
        void gas_molar_density_b(int degree, double y, JetMatrix &Rgb_jet);

        void oil_molar_density(int degree, double y, JetMatrix &Ro_jet);
        void oil_molar_density_b(int degree, const JetMatrix &Ro_jet, double y, JetMatrix &Rob_jet);
        void oil_molar_density_d(int degree, const JetMatrix &Ro_jet, double y, JetMatrix &Rod_jet);
};

#endif // _DEADVOLATILEVOLATILEGASTHERMODYNAMICS_


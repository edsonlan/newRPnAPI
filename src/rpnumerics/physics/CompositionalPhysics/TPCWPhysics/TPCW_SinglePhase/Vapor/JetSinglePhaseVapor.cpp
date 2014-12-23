#include "JetSinglePhaseVapor.h"

JetSinglePhaseVapor::JetSinglePhaseVapor(double mc, double mw, double P) : MC(mc), MW(mw){
    rho_md = new MolarDensity(MOLAR_DENSITY_VAPOR, P);
}

JetSinglePhaseVapor::~JetSinglePhaseVapor(){
    delete rho_md;
}

// Assume JetMatrix has dimension 2 (size = 14)
//
int JetSinglePhaseVapor::rhosigmac_jet(const double yw, const double T, int degree, JetMatrix &m){
    JetMatrix r(2);
    rho_md->rho_jet(yw, T, degree, r);

    if (degree >= 0){
        double rho = r.get(0);
        double rhosigmac = (1.0 - yw)*rho*MC;

        m.set(0, rhosigmac);

        if (degree >= 1){
            double drho_dyw = r.get(0, 0);
            double drho_dT  = r.get(0, 1);

            double drhosigmac_dyw, drhosigmac_dT;

            drhosigmac_dyw = (-rho + (1.0 - yw)*drho_dyw)*MC;

            drhosigmac_dT  = (1.0 - yw)*drho_dT*MC;

            m.set(0, 0, drhosigmac_dyw);
            m.set(0, 1, drhosigmac_dT);

            if (degree == 2){
                double d2rho_dyw2  = r.get(0, 0, 0);
                double d2rho_dywdT = r.get(0, 0, 1);
                double d2rho_dTdyw = r.get(0, 1, 0);
                double d2rho_dT2   = r.get(0, 1, 1);

                double d2rhosigmac_dyw2, d2rhosigmac_dywdT, d2rhosigmac_dTdyw, d2rhosigmac_dT2;

                d2rhosigmac_dyw2  = (-2.0*drho_dyw + (1.0 - yw)*d2rho_dyw2)*MC;
                d2rhosigmac_dywdT = (-drho_dT + (1.0 - yw)*d2rho_dywdT)*MC;
                d2rhosigmac_dTdyw = d2rhosigmac_dywdT;
                d2rhosigmac_dT2   = (1.0 - yw)*d2rho_dT2*MC;

                m.set(0, 0, 0, d2rhosigmac_dyw2);
                m.set(0, 0, 1, d2rhosigmac_dywdT);
                m.set(0, 1, 0, d2rhosigmac_dTdyw);
                m.set(0, 1, 1, d2rhosigmac_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE;
    }
    else return -1; // ABORTED_PROCEDURE
}

int JetSinglePhaseVapor::rhosigmaw_jet(const double yw, const double T, int degree, JetMatrix &m){
    JetMatrix r(2);
    rho_md->rho_jet(yw, T, degree, r);

    if (degree >= 0){
        double rho = r.get(0);
        double rhosigmaw = yw*rho*MW;

        m.set(0, rhosigmaw);
        if (degree >= 1){
            double drho_dyw = r.get(0, 0);
            double drho_dT  = r.get(0, 1);

            double drhosigmaw_dyw, drhosigmaw_dT;

            drhosigmaw_dyw = (rho + yw*drho_dyw)*MW;

            drhosigmaw_dT  = yw*drho_dT*MW;

            m.set(0, 0, drhosigmaw_dyw);
            m.set(0, 1, drhosigmaw_dT);
            if (degree == 2){
                double d2rho_dyw2  = r.get(0, 0, 0);
                double d2rho_dywdT = r.get(0, 0, 1);
                double d2rho_dTdyw = r.get(0, 1, 0);
                double d2rho_dT2   = r.get(0, 1, 1);

                double d2rhosigmaw_dyw2, d2rhosigmaw_dywdT, d2rhosigmaw_dTdyw, d2rhosigmaw_dT2;

                d2rhosigmaw_dyw2  = (2.0*drho_dyw + yw*d2rho_dyw2)*MW;
                d2rhosigmaw_dywdT = (drho_dT + yw*d2rho_dywdT)*MW;
                d2rhosigmaw_dTdyw = d2rhosigmaw_dywdT;
                d2rhosigmaw_dT2   = yw*d2rho_dT2*MW;

                m.set(0, 0, 0, d2rhosigmaw_dyw2);
                m.set(0, 0, 1, d2rhosigmaw_dywdT);
                m.set(0, 1, 0, d2rhosigmaw_dTdyw);
                m.set(0, 1, 1, d2rhosigmaw_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE
    }
    else return -1; // ABORTED_PROCEDURE
}


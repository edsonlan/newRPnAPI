#include "JetSinglePhaseLiquid.h"

JetSinglePhaseLiquid::JetSinglePhaseLiquid(double mc, double mw, double P) : MC(mc), MW(mw){
    rho_md = new MolarDensity(MOLAR_DENSITY_LIQUID, P);
}

JetSinglePhaseLiquid::~JetSinglePhaseLiquid(){
    delete rho_md;
}

// Assume JetMatrix has dimension 2 (size = 14)
//
int JetSinglePhaseLiquid::rhoac_jet(const double xc, const double T, int degree, JetMatrix &m){
    JetMatrix r(2);
    rho_md->rho_jet(xc, T, degree, r);

    if (degree >= 0){
        double rho = r.get(0);
        double rhoac = xc*rho*MC;

        m.set(0, rhoac);

        if (degree >= 1){
            double drho_dxc = r.get(0, 0);
            double drho_dT  = r.get(0, 1);

            double drhoac_dxc, drhoac_dT;

            drhoac_dxc = (rho + xc*drho_dxc)*MC;

            drhoac_dT  = xc*drho_dT*MC;

            m.set(0, 0, drhoac_dxc);
            m.set(0, 1, drhoac_dT);

            if (degree == 2){
                double d2rho_dxc2  = r.get(0, 0, 0);
                double d2rho_dxcdT = r.get(0, 0, 1);
                double d2rho_dTdxc = r.get(0, 1, 0);
                double d2rho_dT2   = r.get(0, 1, 1);

                double d2rhoac_dxc2, d2rhoac_dxcdT, d2rhoac_dTdxc, d2rhoac_dT2;

                d2rhoac_dxc2  = (2.0*drho_dxc + xc*d2rho_dxc2)*MC;
                d2rhoac_dxcdT = (drho_dT + xc*d2rho_dxcdT)*MC;
                d2rhoac_dTdxc = d2rhoac_dxcdT;
                d2rhoac_dT2   = xc*d2rho_dT2*MC;

                m.set(0, 0, 0, d2rhoac_dxc2);
                m.set(0, 0, 1, d2rhoac_dxcdT);
                m.set(0, 1, 0, d2rhoac_dTdxc);
                m.set(0, 1, 1, d2rhoac_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE;
    }
    else return -1; // ABORTED_PROCEDURE
}

int JetSinglePhaseLiquid::rhoaw_jet(const double xc, const double T, int degree, JetMatrix &m){
    JetMatrix r(2); //printf("Before rho_md->rho_jet\n");
    rho_md->rho_jet(xc, T, degree, r);//printf("After rho_md->rho_jet\n");

    if (degree >= 0){
        double rho = r.get(0);
        double rhoaw = (1.0 - xc)*rho*MW;

        m.set(0, rhoaw);
        if (degree >= 1){
            double drho_dxc = r.get(0, 0);
            double drho_dT  = r.get(0, 1);

            double drhoaw_dxc, drhoaw_dT;

            drhoaw_dxc = (-rho + (1.0 - xc)*drho_dxc)*MW;

            drhoaw_dT  = (1.0 - xc)*drho_dT*MW;

            m.set(0, 0, drhoaw_dxc);
            m.set(0, 1, drhoaw_dT);
            if (degree == 2){
                double d2rho_dxc2  = r.get(0, 0, 0);
                double d2rho_dxcdT = r.get(0, 0, 1);
                double d2rho_dTdxc = r.get(0, 1, 0);
                double d2rho_dT2   = r.get(0, 1, 1);

                double d2rhoaw_dxc2, d2rhoaw_dxcdT, d2rhoaw_dTdxc, d2rhoaw_dT2;

                d2rhoaw_dxc2  = (-2.0*drho_dxc + (1.0 - xc)*d2rho_dxc2)*MW;
                d2rhoaw_dxcdT = (-drho_dT + (1.0 - xc)*d2rho_dxcdT)*MW;
                d2rhoaw_dTdxc = d2rhoaw_dxcdT;
                d2rhoaw_dT2   = (1.0 - xc)*d2rho_dT2*MW;

                m.set(0, 0, 0, d2rhoaw_dxc2);
                m.set(0, 0, 1, d2rhoaw_dxcdT);
                m.set(0, 1, 0, d2rhoaw_dTdxc);
                m.set(0, 1, 1, d2rhoaw_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE
    }
    else return -1; // ABORTED_PROCEDURE
}


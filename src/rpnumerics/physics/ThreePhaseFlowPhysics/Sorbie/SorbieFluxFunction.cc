#include "SorbieFluxFunction.h"

SorbieFluxFunction::SorbieFluxFunction(Parameter *muw, Parameter *muo, Parameter *mug, 
                                       Parameter *grw, Parameter *gro, Parameter *grg,
                                       Parameter *vel,
                                       SorbiePermeability *p): FluxFunction(){
    muw_parameter = muw;
    muo_parameter = muo;
    mug_parameter = mug;

    grw_parameter = grw;
    gro_parameter = gro;
    grg_parameter = grg;

    vel_parameter = vel;

    perm = p;
}

SorbieFluxFunction::~SorbieFluxFunction(){
}

int SorbieFluxFunction::jet(const WaveState &w, JetMatrix &m, int degree) const {
    m.resize(2);

    RealVector state(3);
    state(0) = w(0);
    state(1) = w(1);
    state(2) = 1.0 - w(0) - w(1);

    JetMatrix water;
    perm->PermeabilityWater_jet(state, degree, water);

    JetMatrix gas;
    perm->PermeabilityGas_jet(state, degree, gas);

    JetMatrix oil;
    perm->PermeabilityOil_jet(state, degree, oil);

    if (degree >= 0) {
        double sw = w(0);
        double so = w(1);
        double sg = 1.0 - sw - so;

        double muw = muw_parameter->value();
        double muo = muo_parameter->value();
        double mug = mug_parameter->value();

        double grw = grw_parameter->value();
        double gro = gro_parameter->value();
        double grg = grg_parameter->value();

        double vel = vel_parameter->value();

        // I will use frw, grg, gro for Type I, II and III of Sorbie.
        //
        double T = grw + grg + gro;
        double T1 = grw/T;
        double T2 = grg/T;
        double T3 = gro/T;

        double kw = water.get(0);
        double ko = oil.get(0);
        double kg = gas.get(0);

        double lkw = kw / muw; // Water mobility
        double lko = ko / muo; // Oil mobility
        double lkg = kg / mug; // Gas mobility

        double lko1 = (1 - kw - kg) / muo; // Oil mobility on Type I
        double lkg2 = (1 - kw - ko) / mug; // Gas mobility on Type II
        double lkw3 = (1 - kg - ko) / muw; // Water mobility on Type III

        double lk1 = lkw + lkg + lko1;
        double lk2 = lkw + lkg2 + lko;
        double lk3 = lkw3 + lkg + lko;

        double fw = T1*(lkw/lk1) + T2*(lkw/lk2) + T3*(lkw3/lk3);
        double fo = T1*(lko1/lk1) + T2*(lko/lk2) + T3*(lko/lk3);

        m.set(0, fw); // fw
        m.set(1, fo); // fo

        if (degree >= 1) {
            double dkw_dsw = water.get(0, 0);
            double dko_dsw = oil.get(0, 0);
            double dkg_dsw = gas.get(0, 0);

            double ldkw_dsw = dkw_dsw / muw;
            double ldko_dsw = dko_dsw / muo;
            double ldkg_dsw = dkg_dsw / mug;

            double dkw_dso = water.get(0, 1);
            double dko_dso = oil.get(0, 1);
            double dkg_dso = gas.get(0, 1);

            double ldkw_dso = dkw_dso / muw;
            double ldko_dso = dko_dso / muo;
            double ldkg_dso = dkg_dso / mug;

            double ldko1_dsw = -dkw_dsw / muo;
            double ldko1_dso = 0.0;

            double ldkg2_dsw = -dkw_dsw / mug;
            double ldkg2_dso = -dko_dso / mug;

            double ldkw3_dso = -dko_dso / muw;
            double ldkw3_dsw = 0.0;

            double ldk1_dsw = ldkw_dsw + ldkg_dsw + ldko1_dsw;
            double ldk2_dsw = ldkw_dsw + ldkg2_dsw + ldko_dsw;
            double ldk3_dsw = ldkw3_dsw + ldkg_dsw + ldko_dsw;

            double ldk1_dso = ldkw_dso + ldkg_dso + ldko1_dso;
            double ldk2_dso = ldkw_dso + ldkg2_dso + ldko_dso;
            double ldk3_dso = ldkw3_dso + ldkg_dso + ldko_dso;

            double dfw_dsw = T1*( ldkw_dsw - lkw*ldk1_dsw/lk1 ) / lk1
                           + T2*( ldkw_dsw - lkw*ldk2_dsw/lk2 ) / lk2
                           + T3*( ldkw3_dsw - lkw3*ldk3_dsw/lk3 ) / lk3;
            double dfw_dso = T1*( ldkw_dso - lkw*ldk1_dso/lk1 ) / lk1
                           + T2*( ldkw_dso - lkw*ldk2_dso/lk2 ) / lk2
                           + T3*( ldkw3_dso - lkw3*ldk3_dso/lk3 ) / lk3;

            double dfo_dsw = T1*( ldko1_dsw - lko1*ldk1_dsw/lk1 ) / lk1
                           + T2*( ldko_dsw - lko*ldk2_dsw/lk2 ) / lk2
                           + T3*( ldko_dsw - lko*ldk3_dsw/lk3 ) / lk3;
            double dfo_dso = T1*( ldko1_dso - lko1*ldk1_dso/lk1 ) / lk1
                           + T2*( ldko_dso - lko*ldk2_dso/lk2 ) / lk2
                           + T3*( ldko_dso - lko*ldk3_dso/lk3 ) / lk3;


            m.set(0, 0, dfw_dsw);
            m.set(0, 1, dfw_dso);
            m.set(1, 0, dfo_dsw);
            m.set(1, 1, dfo_dso);

            if (degree == 2) {
                // Water.
                //
                double d2kw_dsw2  = water.get(0, 0, 0);
                double d2kw_dswso = water.get(0, 0, 1);
                double d2kw_dso2  = water.get(0, 0, 1);

                double ld2kw_dsw2  = d2kw_dsw2 / muw;
                double ld2kw_dswso = d2kw_dswso / muw;
                double ld2kw_dso2  = d2kw_dso2 / muw;

                // Oil.
                //
                double d2ko_dsw2  = oil.get(0, 0, 0);
                double d2ko_dswso = oil.get(0, 0, 1);
                double d2ko_dso2  = oil.get(0, 0, 1);
            
                double ld2ko_dsw2  = d2ko_dsw2 / muo;
                double ld2ko_dswso = d2ko_dswso / muo;
                double ld2ko_dso2  = d2ko_dso2 / muo;

                // Gas.
                //
                double d2kg_dsw2  = gas.get(0, 0, 0);
                double d2kg_dswso = gas.get(0, 0, 1);
                double d2kg_dso2  = gas.get(0, 0, 1);

                double ld2kg_dsw2  = d2kg_dsw2 / mug;
                double ld2kg_dswso = d2kg_dswso / mug;
                double ld2kg_dso2  = d2kg_dso2 / mug;


                double ld2ko1_dsw2  = -d2kw_dsw2 / muo;
                double ld2ko1_dswso = 0.0;
                double ld2ko1_dso2  = 0.0;

                double ld2kg2_dsw2  = -d2kw_dsw2 / mug;
                double ld2kg2_dswso = -d2kw_dswso / mug;
                double ld2kg2_dso2  = -d2ko_dso2 / mug;

                double ld2kw3_dsw2  = 0.0;
                double ld2kw3_dswso = 0.0;
                double ld2kw3_dso2  = -d2ko_dso2 / muw;

                double ld2k1_dsw2 = ld2kw_dsw2 + ld2kg_dsw2 + ld2ko1_dsw2;
                double ld2k2_dsw2 = ld2kw_dsw2 + ld2kg2_dsw2 + ld2ko_dsw2;
                double ld2k3_dsw2 = ld2kw3_dsw2 + ld2kg_dsw2 + ld2ko_dsw2;

                double ld2k1_dswso = ld2kw_dswso + ld2kg_dswso + ld2ko1_dswso;
                double ld2k2_dswso = ld2kw_dswso + ld2kg2_dswso + ld2ko_dswso;
                double ld2k3_dswso = ld2kw3_dswso + ld2kg_dswso + ld2ko_dswso;

                double ld2k1_dso2 = ldkw_dso + ldkg_dso + ldko1_dso;
                double ld2k2_dso2 = ldkw_dso + ldkg2_dso + ldko_dso;
                double ld2k3_dso2 = ldkw3_dso + ldkg_dso + ldko_dso;

                double d2fw_dsw2  = T1 * ( ld2kw_dsw2 * lk1 - ld2k1_dsw2*lkw - 2.0*ldk1_dsw/lk1*(ldkw_dsw*lk1 - ldk1_dsw*lkw) ) / (lk1*lk1)
                                  + T2 * ( ld2kw_dsw2 * lk2 - ld2k2_dsw2*lkw - 2.0*ldk2_dsw/lk2*(ldkw_dsw*lk2 - ldk2_dsw*lkw) ) / (lk2*lk2)
                                  + T3 * ( ld2kw3_dsw2 * lk3 - ld2k3_dsw2*lkw3 - 2.0*ldk3_dsw/lk2*(ldkw3_dsw*lk3 - ldk3_dsw*lkw3) ) / (lk3*lk3);
                double d2fw_dswso = T1 * ( ld2kw_dswso * lk1 - ld2k1_dswso*lkw - ldk1_dsw/lk1*(ldkw_dso*lk1 - ldk1_dso*lkw) - ldk1_dso/lk1*(ldkw_dsw*lk1 - ldk1_dsw*lkw) ) / (lk1*lk1)
                                  + T2 * ( ld2kw_dswso * lk2 - ld2k2_dswso*lkw - ldk2_dsw/lk2*(ldkw_dso*lk2 - ldk2_dso*lkw) - ldk2_dso/lk2*(ldkw_dsw*lk2 - ldk2_dsw*lkw) ) / (lk2*lk2)
                                  + T3 * ( ld2kw3_dswso * lk3 - ld2k3_dswso*lkw3 - ldk3_dsw/lk2*(ldkw3_dso*lk3 - ldk3_dso*lkw3) - ldk3_dso/lk2*(ldkw3_dsw*lk3 - ldk3_dsw*lkw3) ) / (lk3*lk3);
                double d2fw_dso2  = T1 * ( ld2kw_dso2 * lk1 - ld2k1_dso2*lkw - 2.0*ldk1_dso/lk1*(ldkw_dso*lk1 - ldk1_dso*lkw) ) / (lk1*lk1)
                                  + T2 * ( ld2kw_dso2 * lk2 - ld2k2_dso2*lkw - 2.0*ldk2_dso/lk2*(ldkw_dso*lk2 - ldk2_dso*lkw) ) / (lk2*lk2)
                                  + T3 * ( ld2kw3_dso2 * lk3 - ld2k3_dso2*lkw3 - 2.0*ldk3_dso/lk2*(ldkw3_dso*lk3 - ldk3_dso*lkw3) ) / (lk3*lk3);

                double d2fo_dsw2  = T1 * ( ld2ko1_dsw2 * lk1 - ld2k1_dsw2*lko1 - 2.0*ldk1_dsw/lk1*(ldko1_dsw*lk1 - ldk1_dsw*lko1) ) / (lk1*lk1)
                                  + T2 * ( ld2ko_dsw2 * lk2 - ld2k2_dsw2*lko - 2.0*ldk2_dsw/lk2*(ldko_dsw*lk2 - ldk2_dsw*lko) ) / (lk2*lk2)
                                  + T3 * ( ld2ko_dsw2 * lk3 - ld2k3_dsw2*lko - 2.0*ldk3_dsw/lk2*(ldko_dsw*lk3 - ldk3_dsw*lko) ) / (lk3*lk3);
                double d2fo_dswso = T1 * ( ld2ko1_dswso * lk1 - ld2k1_dswso*lko1 - ldk1_dsw/lk1*(ldko1_dso*lk1 - ldk1_dso*lko1) - ldk1_dso/lk1*(ldko1_dsw*lk1 - ldk1_dsw*lko1) ) / (lk1*lk1)
                                  + T2 * ( ld2ko_dswso * lk2 - ld2k2_dswso*lko - ldk2_dsw/lk2*(ldko_dso*lk2 - ldk2_dso*lko) - ldk2_dso/lk2*(ldko_dsw*lk2 - ldk2_dsw*lko) ) / (lk2*lk2)
                                  + T3 * ( ld2ko_dswso * lk3 - ld2k3_dswso*lko - ldk3_dsw/lk2*(ldko_dso*lk3 - ldk3_dso*lko) - ldk3_dso/lk2*(ldko_dsw*lk3 - ldk3_dsw*lko) ) / (lk3*lk3);
                double d2fo_dso2  = T1 * ( ld2ko1_dso2 * lk1 - ld2k1_dso2*lko1 - 2.0*ldk1_dso/lk1*(ldko1_dso*lk1 - ldk1_dso*lko1) ) / (lk1*lk1)
                                  + T2 * ( ld2ko_dso2 * lk2 - ld2k2_dso2*lko - 2.0*ldk2_dso/lk2*(ldko_dso*lk2 - ldk2_dso*lko) ) / (lk2*lk2)
                                  + T3 * ( ld2ko_dso2 * lk3 - ld2k3_dso2*lko - 2.0*ldk3_dso/lk2*(ldko_dso*lk3 - ldk3_dso*lko) ) / (lk3*lk3);

                m.set(0, 0, 0, d2fw_dsw2);
                m.set(0, 0, 1, d2fw_dswso);
                m.set(0, 1, 0, d2fw_dswso);
                m.set(0, 1, 1, d2fw_dso2);

                m.set(1, 0, 0, d2fo_dsw2);
                m.set(1, 0, 1, d2fo_dswso);
                m.set(1, 1, 0, d2fo_dswso);
                m.set(1, 1, 1, d2fo_dso2);
            }
        }
    }

    return 2; //SUCCESSFUL_PROCEDURE;
}


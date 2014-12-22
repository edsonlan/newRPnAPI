#include "Flux2Comp2PhasesAdimensionalized.h"

Flux2Comp2PhasesAdimensionalized::Flux2Comp2PhasesAdimensionalized(Parameter *abs_perm, Parameter *sin_beta, 
                                                                   Parameter *cnw, Parameter *cng,
                                                                   Parameter *expw, Parameter *expg,
                                                                   bool has_grav, bool has_hor,
                                                                   Thermodynamics *td): FluxFunction()
                                                                   #ifdef JETTESTER_ENABLED_TPCWACCUMULATION
                                                                   , TestableJet()
                                                                   #endif  
{
    abs_perm_parameter = abs_perm;
    sin_beta_parameter = sin_beta;
    cnw_parameter      = cnw;
    cng_parameter      = cng;
    expw_parameter     = expw;
    expg_parameter     = expg;

    has_gravity    = has_grav;
    has_horizontal = has_hor;

    const_gravity = 9.8;

    TD = td;

    FH = new FracFlow2PhasesHorizontalAdimensionalized(this);
    FV = new FracFlow2PhasesVerticalAdimensionalized(this);
    reducedFlux = new ReducedFlux2Comp2PhasesAdimensionalized(this);
}

Flux2Comp2PhasesAdimensionalized::~Flux2Comp2PhasesAdimensionalized() {
    delete reducedFlux;
    delete FV;
    delete FH;
}

//void Flux2Comp2PhasesAdimensionalized::fluxParams(const FluxParams & param) {

//    FluxFunction::fluxParams(param);




//    if (param.component(2) == 1.0) {
//        has_gravity = true;
//    } else {
//        has_gravity = false;
//    }

//    if (param.component(3) == 1.0) {
//        has_horizontal = true;
//    } else {
//        has_horizontal = false;
//    }

//    abs_perm = param.component(0);
//    sin_beta = param.component(1);


//    const_gravity = 9.8;

//    cnw = param.component(4);
//    cng = param.component(5);
//    expw = param.component(6);
//    expg = param.component(7);

//    grav = abs_perm * sin_beta*const_gravity;



//}

int Flux2Comp2PhasesAdimensionalized::ReducedFlux2Comp2PhasesAdimensionalized::jet(const WaveState &w, JetMatrix &m, int degree) const {
    double s = w(0); // s_{sigma} = sg in FracFlow2PhasesHorizontal & FracFlow2PhasesVertical
    double Theta = w(1);

    // Some auxiliary variables
    JetMatrix Hrj(1);
    fluxComplete_->TD->RockEnthalpyVol_jet(Theta, degree, Hrj);

    JetMatrix Haj(1);
    fluxComplete_->TD->AqueousEnthalpyVol_jet(Theta, degree, Haj);

    JetMatrix Hsij(1);
    fluxComplete_->TD->SuperCriticEnthalpyVol_jet(Theta, degree, Hsij);

    JetMatrix rhosicj(1);
    fluxComplete_->TD->Rhosic_jet(Theta, degree, rhosicj);

    JetMatrix rhosiwj(1);
    fluxComplete_->TD->Rhosiw_jet(Theta, degree, rhosiwj);

    JetMatrix rhoacj(1);
    fluxComplete_->TD->Rhoac_jet(Theta, degree, rhoacj);

    JetMatrix rhoawj(1);
    fluxComplete_->TD->Rhoaw_jet(Theta, degree, rhoawj);

    // Fractional Flow Function.
    JetMatrix horizontal(2);
    fluxComplete_->FH->Diff_FracFlow2PhasesHorizontalAdimensionalized(1. - s, Theta, degree, horizontal);

    if (degree >= 0) {
        double Hr     = Hrj.get(0);
        double Ha     = Haj.get(0);
        double Hsi    = Hsij.get(0);
        double rhosic = rhosicj.get(0);
        double rhosiw = rhosiwj.get(0);
        double rhoac  = rhoacj.get(0);
        double rhoaw  = rhoawj.get(0);

        double f = horizontal.get(0);

        double out0 = (rhosic * f + rhoac * (1.0 - f));
        double out1 = (rhosiw * f + rhoaw * (1.0 - f));
        double out2 = (Hsi * f + Ha * (1.0 - f));

        m.set(0, out0);
        m.set(1, out1);
        m.set(2, out2);

        if (degree >= 1) {
            double d_Hr     = Hrj.get(0, 0);
            double d_Ha     = Haj.get(0, 0);
            double d_Hsi    = Hsij.get(0, 0);
            double d_rhosic = rhosicj.get(0, 0);
            double d_rhosiw = rhosiwj.get(0, 0);
            double d_rhoac  = rhoacj.get(0, 0);
            double d_rhoaw  = rhoawj.get(0, 0);

            double df_ds     = horizontal.get(0, 0);
            double df_dTheta = horizontal.get(0, 1);

            double out00 = (rhosic - rhoac) * df_ds;
            double out01 = (d_rhosic * f + d_rhoac * (1.0 - f) + (rhosic - rhoac) * df_dTheta);

            double out10 = (rhosiw - rhoaw) * df_ds;
            double out11 = (d_rhosiw * f + d_rhoaw * (1.0 - f) + (rhosiw - rhoaw) * df_dTheta);

            double out20 = (Hsi - Ha) * df_ds;
            double out21 = (d_Hsi * f + d_Ha * (1.0 - f) + (Hsi - Ha) * df_dTheta);

            m.set(0, 0, out00);
            m.set(0, 1, out01);

            m.set(1, 0, out10);
            m.set(1, 1, out11);

            m.set(2, 0, out20);
            m.set(2, 1, out21);

            if (degree >= 2) {
                double d2_Hr     = Hrj.get(0, 0, 0);
                double d2_Ha     = Haj.get(0, 0, 0);
                double d2_Hsi    = Hsij.get(0, 0, 0);
                double d2_rhosic = rhosicj.get(0, 0, 0);
                double d2_rhosiw = rhosiwj.get(0, 0, 0);
                double d2_rhoac  = rhoacj.get(0, 0, 0);
                double d2_rhoaw  = rhoawj.get(0, 0, 0);

                double d2f_ds2      = horizontal.get(0, 0, 0);
                double d2f_dsdTheta = horizontal.get(0, 0, 1);
                double d2f_dThetads = horizontal.get(0, 1, 0);
                double d2f_dTheta2  = horizontal.get(0, 1, 1);

                double out000 = (rhosic - rhoac) * d2f_ds2;
                double out001 = ((d_rhosic - d_rhoac) * df_ds + (rhosic - rhoac) * d2f_dsdTheta);

                double out010 = out001; // Mixed partial
                double out011 = (d2_rhosic * f + d2_rhoac * (1.0 - f) + 2 * (d_rhosic - d_rhoac) * df_dTheta + (rhosic - rhoac) * d2f_dTheta2);

                double out100 = (rhosiw - rhoaw) * d2f_ds2;
                double out101 = ((d_rhosiw - d_rhoaw) * df_ds + (rhosiw - rhoaw) * d2f_dsdTheta);

                double out110 = out101; // Mixed partial
                double out111 = (d2_rhosiw * f + d2_rhoaw * (1.0 - f) + 2 * (d_rhosiw - d_rhoaw) * df_dTheta + (rhosiw - rhoaw) * d2f_dTheta2);

                double out200 = (Hsi - Ha) * d2f_ds2;
                double out201 = ((d_Hsi - d_Ha) * df_ds + (Hsi - Ha) * d2f_dsdTheta);

                double out210 = out201; // Mixed partial
                double out211 = (d2_Hsi * f + d2_Ha * (1.0 - f) + 2 * (d_Hsi - d_Ha) * df_dTheta + (Hsi - Ha) * d2f_dTheta2);

                m.set(0, 0, 0, out000);
                m.set(0, 0, 1, out001);
                m.set(0, 1, 0, out010);
                m.set(0, 1, 1, out011);

                m.set(1, 0, 0, out100);
                m.set(1, 0, 1, out101);
                m.set(1, 1, 0, out110);
                m.set(1, 1, 1, out111);

                m.set(2, 0, 0, out200);
                m.set(2, 0, 1, out201);
                m.set(2, 1, 0, out210);
                m.set(2, 1, 1, out211);

            }
        }
    }

    return 2; //SUCCESSFUL_PROCEDURE;
}

Flux2Comp2PhasesAdimensionalized::ReducedFlux2Comp2PhasesAdimensionalized::ReducedFlux2Comp2PhasesAdimensionalized(Flux2Comp2PhasesAdimensionalized * outer){
    fluxComplete_ = outer;
}

Flux2Comp2PhasesAdimensionalized::FracFlow2PhasesHorizontalAdimensionalized::FracFlow2PhasesHorizontalAdimensionalized(Flux2Comp2PhasesAdimensionalized * outer) {
    fluxComplete_ = outer;
}

Flux2Comp2PhasesAdimensionalized::FracFlow2PhasesVerticalAdimensionalized::FracFlow2PhasesVerticalAdimensionalized(Flux2Comp2PhasesAdimensionalized * outer) {

    fluxComplete_ = outer;

}

//RpFunction * Flux2Comp2PhasesAdimensionalized::ReducedFlux2Comp2PhasesAdimensionalized::clone() const {

//}

int Flux2Comp2PhasesAdimensionalized::jet(const WaveState &w, JetMatrix &m, int degree) const {
    double s     = w(0); // s_{sigma} = sg in FracFlow2PhasesHorizontal & FracFlow2PhasesVertical
    double Theta = w(1);
    double U     = w(2);
    
    double abs_perm = abs_perm_parameter->value();
    double sin_beta = sin_beta_parameter->value();
    double cnw = cnw_parameter->value();
    double cng = cng_parameter->value();
    double expw = expw_parameter->value();
    double expg = expg_parameter->value();

    double grav = abs_perm*sin_beta*const_gravity;


    m.resize(3);    

    JetMatrix Hrj(1);
    TD->RockEnthalpyVol_jet(Theta, degree, Hrj);

    JetMatrix Haj(1);
    TD->AqueousEnthalpyVol_jet(Theta, degree, Haj);

    JetMatrix Hsij(1);
    TD->SuperCriticEnthalpyVol_jet(Theta, degree, Hsij);

    JetMatrix rhosicj(1);
    TD->Rhosic_jet(Theta, degree, rhosicj);

    JetMatrix rhosiwj(1);
    TD->Rhosiw_jet(Theta, degree, rhosiwj);

    JetMatrix rhoacj(1);
    TD->Rhoac_jet(Theta, degree, rhoacj);

    JetMatrix rhoawj(1);
    TD->Rhoaw_jet(Theta, degree, rhoawj);

    // Fill here, since below will be more complicated.
    double Hr     = Hrj.get(0);
    double Ha     = Haj.get(0);
    double Hsi    = Hsij.get(0);
    double rhosic = rhosicj.get(0);
    double rhosiw = rhosiwj.get(0);
    double rhoac  = rhoacj.get(0);
    double rhoaw  = rhoawj.get(0);

    double d_Hr     = Hrj.get(0, 0);
    double d_Ha     = Haj.get(0, 0);
    double d_Hsi    = Hsij.get(0, 0);
    double d_rhosic = rhosicj.get(0, 0);
    double d_rhosiw = rhosiwj.get(0, 0);
    double d_rhoac  = rhoacj.get(0, 0);
    double d_rhoaw  = rhoawj.get(0, 0);

    double d2_Hr     = Hrj.get(0, 0, 0);
    double d2_Ha     = Haj.get(0, 0, 0);
    double d2_Hsi    = Hsij.get(0, 0, 0);
    double d2_rhosic = rhosicj.get(0, 0, 0);
    double d2_rhosiw = rhosiwj.get(0, 0, 0);
    double d2_rhoac  = rhoacj.get(0, 0, 0);
    double d2_rhoaw  = rhoawj.get(0, 0, 0);

    // Output
    double out0 = 0.0, out1 = 0.0, out2 = 0.0;

    double out00 = 0.0, out01 = 0.0, out02 = 0.0;
    double out10 = 0.0, out11 = 0.0, out12 = 0.0;
    double out20 = 0.0, out21 = 0.0, out22 = 0.0;

    double out000 = 0.0, out001 = 0.0, out002 = 0.0;
    double out010 = 0.0, out011 = 0.0, out012 = 0.0;
    double out020 = 0.0, out021 = 0.0, out022 = 0.0;

    double out100 = 0.0, out101 = 0.0, out102 = 0.0;
    double out110 = 0.0, out111 = 0.0, out112 = 0.0;
    double out120 = 0.0, out121 = 0.0, out122 = 0.0;

    double out200 = 0.0, out201 = 0.0, out202 = 0.0;
    double out210 = 0.0, out211 = 0.0, out212 = 0.0;
    double out220 = 0.0, out221 = 0.0, out222 = 0.0;

    // Begin of pure horizontal

    if (has_horizontal) {

        double f, df_ds, df_dTheta, d2f_ds2, d2f_dsdTheta, d2f_dThetads, d2f_dTheta2; // f=f_{sigma}, s=s_{sigma}
        JetMatrix horizontal(2);

        FH->Diff_FracFlow2PhasesHorizontalAdimensionalized(1. - s, Theta, degree, horizontal);

        if (degree >= 0) {
            f = horizontal.get(0);
            if (degree >= 1) {
                df_ds = horizontal.get(0, 0);
                df_dTheta = horizontal.get(0, 1);
                if (degree >= 2) {
                    d2f_ds2 = horizontal.get(0, 0, 0);
                    d2f_dsdTheta = horizontal.get(0, 0, 1);
                    d2f_dThetads = horizontal.get(0, 1, 0);
                    d2f_dTheta2 = horizontal.get(0, 1, 1);
                }
            }
        }

        if (degree >= 0) {
            out0 = U * (rhosic * f + rhoac * (1.0 - f));
            out1 = U * (rhosiw * f + rhoaw * (1.0 - f));
            out2 = U * (Hsi * f + Ha * (1.0 - f));

            if (degree >= 1) {
                out00 = U * (rhosic - rhoac) * df_ds;
                out01 = U * (d_rhosic * f + d_rhoac * (1.0 - f) + (rhosic - rhoac) * df_dTheta);
                out02 = (rhosic * f + rhoac * (1.0 - f)); // dF1_dU

                out10 = U * (rhosiw - rhoaw) * df_ds;
                out11 = U * (d_rhosiw * f + d_rhoaw * (1.0 - f) + (rhosiw - rhoaw) * df_dTheta);
                out12 = (rhosiw * f + rhoaw * (1.0 - f)); // dF2_dU

                out20 = U * (Hsi - Ha) * df_ds;
                out21 = U * (d_Hsi * f + d_Ha * (1.0 - f) + (Hsi - Ha) * df_dTheta);
                out22 = (Hsi * f + Ha * (1.0 - f)); // dF3_dU

                if (degree >= 2) {
                    out000 = U * (rhosic - rhoac) * d2f_ds2;
                    out001 = U * ((d_rhosic - d_rhoac) * df_ds + (rhosic - rhoac) * d2f_dsdTheta);
                    out002 = ((rhosic - rhoac) * df_ds); // d2F1_dsdU

                    out010 = out001; // Mixed partial
                    out011 = U * (d2_rhosic * f + d2_rhoac * (1.0 - f) + 2 * (d_rhosic - d_rhoac) * df_dTheta + (rhosic - rhoac) * d2f_dTheta2);
                    out012 = (d_rhosic * f + d_rhoac * (1.0 - f) + (rhosic - rhoac) * df_dTheta); // d2F1_dThetadU

                    out100 = U * (rhosiw - rhoaw) * d2f_ds2;
                    out101 = U * ((d_rhosiw - d_rhoaw) * df_ds + (rhosiw - rhoaw) * d2f_dsdTheta);
                    out102 = (rhosiw - rhoaw) * df_ds; // d2F2_dsdU

                    out110 = out101; // Mixed partial
                    out111 = U * (d2_rhosiw * f + d2_rhoaw * (1.0 - f) + 2 * (d_rhosiw - d_rhoaw) * df_dTheta + (rhosiw - rhoaw) * d2f_dTheta2);
                    out112 = (d_rhosiw * f + d_rhoaw * (1.0 - f) + (rhosiw - rhoaw) * df_dTheta); // d2F2_dThetadU

                    out120 = out102;
                    out121 = out112;
                    out122 = 0.; // d2F2_dU2

                    out200 = U * (Hsi - Ha) * d2f_ds2;
                    out201 = U * ((d_Hsi - d_Ha) * df_ds + (Hsi - Ha) * d2f_dsdTheta);
                    out202 = (Hsi - Ha) * df_ds; // d2F3_dsdU

                    out210 = out201; // Mixed partial
                    out211 = U * (d2_Hsi * f + d2_Ha * (1.0 - f) + 2 * (d_Hsi - d_Ha) * df_dTheta + (Hsi - Ha) * d2f_dTheta2);
                    out212 = d_Hsi * f + d_Ha * (1.0 - f) + (Hsi - Ha) * df_dTheta; // d2F3_dThetadU

                    out220 = (Hsi - Ha) * df_ds;
                    out221 = d_Hsi * f + d_Ha * (1.0 - f) + (Hsi - Ha) * df_dTheta;
                    out222 = 0.; // d2F3_dU2
                }
            }
        }
    } // End of pure horizontal

    // Begin of pure gravity
    if (has_gravity) {
        double Z, dZ_ds, dZ_dTheta, d2Z_ds2, d2Z_dsdTheta, d2Z_dThetads, d2Z_dTheta2;
        JetMatrix vertical(2);

        FV->Diff_FracFlow2PhasesVerticalAdimensionalized(1. - s, Theta, degree, vertical);

        if (degree >= 0) {
            Z = vertical.get(0);
            if (degree >= 1) {
                dZ_ds = vertical.get(0, 0);
                dZ_dTheta = vertical.get(0, 1);
                if (degree >= 2) {
                    d2Z_ds2 = vertical.get(0, 0, 0);
                    d2Z_dsdTheta = vertical.get(0, 0, 1);
                    d2Z_dThetads = vertical.get(0, 1, 0);
                    d2Z_dTheta2 = vertical.get(0, 1, 1);
                }
            }
        }

        double grhoa_rhosi = grav * ((rhoac + rhoaw) - (rhosic + rhosiw));
        double d_grhoa_rhosi = grav * ((d_rhoac + d_rhoaw) - (d_rhosic + d_rhosiw));
        double d2_grhoa_rhosi = grav * ((d2_rhoac + d2_rhoaw) - (d2_rhosic + d2_rhosiw));

        double rhosic_rhoac = rhosic - rhoac;
        double d_rhosic_rhoac = d_rhosic - d_rhoac;
        double d2_rhosic_rhoac = d2_rhosic - d2_rhoac;

        double rhosiw_rhoaw = rhosiw - rhoaw;
        double d_rhosiw_rhoaw = d_rhosiw - d_rhoaw;
        double d2_rhosiw_rhoaw = d2_rhosiw - d2_rhoaw;

        double Hsi_Ha = Hsi - Ha;
        double d_Hsi_Ha = d_Hsi - d_Ha;
        double d2_Hsi_Ha = d2_Hsi - d2_Ha;

        if (degree >= 0) {
            out0 += grhoa_rhosi * rhosic_rhoac*Z; // F1
            out1 += grhoa_rhosi * rhosiw_rhoaw*Z; // F2
            out2 += grhoa_rhosi * Hsi_Ha*Z; // F3

            if (degree >= 1) {
                out00 += grhoa_rhosi * rhosic_rhoac*dZ_ds; // dF1_ds
                out01 += d_grhoa_rhosi * rhosic_rhoac * Z + grhoa_rhosi * d_rhosic_rhoac * Z + grhoa_rhosi * rhosic_rhoac*dZ_dTheta; // dF1_dTheta
                //  out02 += 0; // dF1_du

                out10 += grhoa_rhosi * rhosiw_rhoaw*dZ_ds; // dF2_ds
                out11 += d_grhoa_rhosi * rhosiw_rhoaw * Z + grhoa_rhosi * d_rhosiw_rhoaw * Z + grhoa_rhosi * rhosiw_rhoaw*dZ_dTheta; // dF2_dTheta
                //  out12 += 0; // dF2_du

                out20 += grhoa_rhosi * Hsi_Ha*dZ_ds; // dF3_ds
                out21 += d_grhoa_rhosi * Hsi_Ha * Z + grhoa_rhosi * d_Hsi_Ha * Z + grhoa_rhosi * Hsi_Ha*dZ_dTheta; // dF3_dTheta
                //  out22 += 0; // dF3_du

                if (degree >= 2) {
                    out000 += grhoa_rhosi * rhosic_rhoac*d2Z_ds2; // d2F1_ds2
                    out001 += d_grhoa_rhosi * rhosic_rhoac * dZ_ds +
                            grhoa_rhosi * d_rhosic_rhoac * dZ_ds +
                            grhoa_rhosi * rhosic_rhoac*d2Z_dsdTheta; // d2F1_dsdTheta
                    // out002 += 0; // d2F1_dsdu

                    out010 += d_grhoa_rhosi * rhosic_rhoac * dZ_ds +
                            grhoa_rhosi * d_rhosic_rhoac * dZ_ds +
                            grhoa_rhosi * rhosic_rhoac*d2Z_dsdTheta; // d2F1_dThetads
                    out011 += d2_grhoa_rhosi * rhosic_rhoac * Z +
                            d_grhoa_rhosi * d_rhosic_rhoac * Z +
                            d_grhoa_rhosi * rhosic_rhoac * dZ_dTheta +
                            d_grhoa_rhosi * d_rhosic_rhoac * Z +
                            grhoa_rhosi * d2_rhosic_rhoac * Z +
                            grhoa_rhosi * d_rhosic_rhoac * dZ_dTheta +
                            d_grhoa_rhosi * rhosic_rhoac * dZ_dTheta +
                            grhoa_rhosi * d_rhosic_rhoac * dZ_dTheta +
                            grhoa_rhosi * rhosic_rhoac*d2Z_dTheta2; // d2F1_dTheta2
                    //  out012 += 0; // d2F1_dThetadu

                    //  out020 += 0; // d2F1_duds
                    //  out021 += 0; // d2F1_dudTheta
                    //  out022 += 0; // d2F1_du2

                    out100 += grhoa_rhosi * rhosiw_rhoaw*d2Z_ds2; // d2F2_ds2
                    out101 += d_grhoa_rhosi * rhosiw_rhoaw * dZ_ds +
                            grhoa_rhosi * d_rhosiw_rhoaw * dZ_ds +
                            grhoa_rhosi * rhosiw_rhoaw*d2Z_dsdTheta; // d2F2_dsdTheta
                    //  out102 += 0; // d2F2_dsdu

                    out110 += d_grhoa_rhosi * rhosiw_rhoaw * dZ_ds +
                            grhoa_rhosi * d_rhosiw_rhoaw * dZ_ds +
                            grhoa_rhosi * rhosiw_rhoaw*d2Z_dsdTheta; // d2F2_dThetads
                    out111 += d2_grhoa_rhosi * rhosiw_rhoaw * Z +
                            d_grhoa_rhosi * d_rhosiw_rhoaw * Z +
                            d_grhoa_rhosi * rhosiw_rhoaw * dZ_dTheta +
                            d_grhoa_rhosi * d_rhosiw_rhoaw * Z +
                            grhoa_rhosi * d2_rhosiw_rhoaw * Z +
                            grhoa_rhosi * d_rhosiw_rhoaw * dZ_dTheta +
                            d_grhoa_rhosi * rhosiw_rhoaw * dZ_dTheta +
                            grhoa_rhosi * d_rhosiw_rhoaw * dZ_dTheta +
                            grhoa_rhosi * rhosiw_rhoaw*d2Z_dTheta2; // d2F2_dTheta2
                    //  out112 += 0; // d2F2_dThetadu

                    //  out120 += 0; // d2F2_duds
                    //  out121 += 0; // d2F2_dudTheta
                    //  out122 += 0; // d2F2_du2

                    out200 += grhoa_rhosi * Hsi_Ha*d2Z_ds2; // d2F3_ds2
                    out201 += d_grhoa_rhosi * Hsi_Ha * dZ_ds +
                            grhoa_rhosi * d_Hsi_Ha * dZ_ds +
                            grhoa_rhosi * Hsi_Ha*d2Z_dsdTheta; // d2F3_dsdTheta
                    //  out202 += 0; // d2F3_dsdu

                    out210 += d_grhoa_rhosi * Hsi_Ha * dZ_ds +
                            grhoa_rhosi * d_Hsi_Ha * dZ_ds +
                            grhoa_rhosi * Hsi_Ha*d2Z_dThetads; // d2F3_dThetads
                    out211 += d2_grhoa_rhosi * Hsi_Ha * Z +
                            d_grhoa_rhosi * d_Hsi_Ha * Z +
                            d_grhoa_rhosi * Hsi_Ha * dZ_dTheta +
                            d_grhoa_rhosi * d_Hsi_Ha * Z +
                            grhoa_rhosi * d2_Hsi_Ha * Z +
                            grhoa_rhosi * d_Hsi_Ha * dZ_dTheta +
                            d_grhoa_rhosi * Hsi_Ha * dZ_dTheta +
                            grhoa_rhosi * d_Hsi_Ha * dZ_dTheta +
                            grhoa_rhosi * Hsi_Ha*d2Z_dTheta2; // d2F3_dTheta2
                    //  out212 += 0; // d2F3_dThetadu

                    //  out220 += 0; // d2F3_duds
                    //  out221 += 0; // d2F3_dudTheta
                    //  out222 += 0; // d2F3_du2
                }
            }
        }
    } // End of pure gravity

    if (degree >= 0) {
        m.set(0, out0);
        m.set(1, out1);
        m.set(2, out2);

        if (degree >= 1) {
            m.set(0, 0, out00);
            m.set(0, 1, out01);
            m.set(0, 2, out02);

            m.set(1, 0, out10);
            m.set(1, 1, out11);
            m.set(1, 2, out12);

            m.set(2, 0, out20);
            m.set(2, 1, out21);
            m.set(2, 2, out22);

            //            for (int i = 0; i < 3; i++){
            //                for (int j = 0; j < 3; j++){
            //                    printf("ff(%d, %d) = %g\n", i, j, m(i, j));
            //                }
            //            }

            if (degree >= 2) {
                m.set(0, 0, 0, out000);
                m.set(0, 0, 1, out001);
                m.set(0, 0, 2, out002);
                m.set(0, 1, 0, out010);
                m.set(0, 1, 1, out011);
                m.set(0, 1, 2, out012);
                m.set(0, 2, 0, out020);
                m.set(0, 2, 1, out021);
                m.set(0, 2, 2, out022);

                m.set(1, 0, 0, out100);
                m.set(1, 0, 1, out101);
                m.set(1, 0, 2, out102);
                m.set(1, 1, 0, out110);
                m.set(1, 1, 1, out111);
                m.set(1, 1, 2, out112);
                m.set(1, 2, 0, out120);
                m.set(1, 2, 1, out121);
                m.set(1, 2, 2, out122);

                m.set(2, 0, 0, out200);
                m.set(2, 0, 1, out201);
                m.set(2, 0, 2, out202);
                m.set(2, 1, 0, out210);
                m.set(2, 1, 1, out211);
                m.set(2, 1, 2, out212);
                m.set(2, 2, 0, out220);
                m.set(2, 2, 1, out221);
                m.set(2, 2, 2, out222);
            }
        }
    }

    return 2; //SUCCESSFUL_PROCEDURE;
}

Thermodynamics * Flux2Comp2PhasesAdimensionalized::getThermo() const {
    return TD;
}

Flux2Comp2PhasesAdimensionalized::FracFlow2PhasesHorizontalAdimensionalized * Flux2Comp2PhasesAdimensionalized::getHorizontalFlux()const {
    return FH;
}

Flux2Comp2PhasesAdimensionalized::FracFlow2PhasesVerticalAdimensionalized * Flux2Comp2PhasesAdimensionalized::getVerticalFlux()const {
    return FV;
}

Flux2Comp2PhasesAdimensionalized::ReducedFlux2Comp2PhasesAdimensionalized * Flux2Comp2PhasesAdimensionalized::getReducedFlux()const {
    return reducedFlux;
}

int Flux2Comp2PhasesAdimensionalized::FracFlow2PhasesHorizontalAdimensionalized::Diff_FracFlow2PhasesHorizontalAdimensionalized(double sw, double Theta, int degree, JetMatrix &m) {
    // Esta funcao retorna a funcao de fluxo da agua e seu jato;
    // esta correta para uso do Helmut apenas
    // ja tiramos grw, grg de todos os lugares e as variaveis que as usavam


    double T = fluxComplete_->TD->Theta2T(Theta);

    //    double T = Flux2Comp2PhasesAdimensionalized::getThermo()->Theta2T(Theta);

    // Auxiliary variables
    double sg, law, lag;
    double dlaw_dsw, dlag_dsw;
    double d2law_dsw2, d2lag_dsw2;
    double la, dla_dsw, d2la_dsw2;
    double dflw_dsw, dflw_dT; // Derivatives of the Fractional Flow of water
    double d2flw_dsw2; // d2flw_dsw2 // Second derivative of Fractional Flow of water

    double flw; // Fractional Flow for water

    // Some variables to facilitate the operations with powers
    double pow0, pow1, pow2;

    double nuw, nug; // nuw= 1/muw , nug=1/mug
    double dnuw_dT, d2nuw_dT2; // This values are assigned using TD
    // dnuw_dT = dnuw_dT(T).
    // d2nuw_dT2 = d2nuw_dT2(T).

    double dnug_dT, d2nug_dT2; // This values are assigned using TD
    // dnug_dT = dnug_dT(T),
    // d2nug_dT2 = d2nug_dT2(T).

    double muw, dmuw_dT, d2muw_dT2;
    double mug, dmug_dT, d2mug_dT2;

    double dlaw_dT, d2law_dT2; // Temperature derivatives of law (are not passed as parameters).

    double dlag_dT, d2lag_dT2; // Temperature derivatives of lag (are not passed as parameters).

    double d2law_dswdT; // Mixed derivative (it is not passed as parameter).
    double d2lag_dswdT; // Mixed derivative (it is not passed as parameter).

    double dla_dT;
    double d2flw_dswdT;
    double d2la_dswdT;
    double d2flw_dT2;
    double d2la_dT2;

    // Finally we assign sg
    sg = 1. - sw;

    fluxComplete_->TD->inv_muw(T, nuw, dnuw_dT, d2nuw_dT2);
    fluxComplete_->TD->inv_mug(T, nug, dnug_dT, d2nug_dT2);

    fluxComplete_->TD->muw(T, muw, dmuw_dT, d2muw_dT2);
    fluxComplete_->TD->mug(T, mug, dmug_dT, d2mug_dT2);

    //--------------------------------------------------------------------

    //

    double cnw  = fluxComplete_->cnw_parameter->value();
    double cng  = fluxComplete_->cng_parameter->value();
    double expw = fluxComplete_->expw_parameter->value();
    double expg = fluxComplete_->expg_parameter->value();

    // Before evaluating proper
    if (sw < cnw) {
        law = 0.;
        dlaw_dsw = 0.;
        d2law_dsw2 = 0.;

        dlaw_dT = 0.;
        d2law_dT2 = 0.;

        d2law_dswdT = 0.; // remember that d2law_dswdT is equal to d2law_dTdsw
    } 
    else {
        double temp = (sw - cnw) / (1. - cnw - cng);

        pow2 = 0.5 * pow(temp, expw - 2.); // 0.5*( (sw-s_{connate water})/(1-s_{connate water}-s_{con gas}) )^{expw-2}
        // This model coincides with the formulae in Helmut's Thesis. 0.5 is the re-scaling factor.

        pow1 = temp*pow2; // 0.5*( (sw-s_{connate water})/(1-s_{connate water}-s_{con gas}) )^{expw-1}
        // Observe that pow1 = pow((sw - cnw)/(1.-cnw-cng),expw-1) ;

        pow0 = temp*pow1; // 0.5*( (sw-s_{connate water})/(1-s_{connate water-s_{con gas}}) )^{expw}
        // Observe that pow0 = pow((sw - cnw)/(1.-cnw-cng),expw)


        law = pow0*nuw; // k_{rw}/mu_{w}(T) This is equivalent to k_{ra}/mu_{a} in Helmut's thesis.
        dlaw_dsw = expw * pow1*nuw; // expw*( (sw-s_{connate water })/(1.-cnw-cng) )^{expw-1}/mu_{w}  derivation w.r.t sw
        d2law_dsw2 = (expw - 1.) * expw * pow2*nuw; // (expw - 1.) * expw*( (sw-s_{connate water })/(1.-cnw-cng) )^{expw-2}/mu_{w}

        dlaw_dT = pow0*dnuw_dT; // the value of the derivatives of muw
        // will be passed by the viscosity function.
        d2law_dT2 = pow0*d2nuw_dT2; //

        d2law_dswdT = expw * pow1*dnuw_dT;
    }

    if (sg < cng) {
        lag = 0.;
        dlag_dsw = 0.;
        d2lag_dsw2 = 0.;

        dlag_dT = 0.;
        d2lag_dT2 = 0.;

        d2lag_dswdT = 0.;
    } else {

        double temp = (sg - cng) / (1. - cnw - cng);
        pow2 = 0.95 * pow(temp, expg - 2.); // ( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-2}
        pow1 = temp*pow2; // ( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-1} Observe that pow1 = pow(temp, expg-1)
        pow0 = temp*pow1; // ( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg}   Observe that pow0 = pow(temp, expg) ;

        lag = pow0*nug; //  k_{rg}/mu_{g}(T)  this viscosity is function of temperature   ;
        dlag_dsw = -expg * pow1*nug; //  - expg*( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-1}/mu_{g}
        d2lag_dsw2 = (expg - 1.) * expg * pow2*nug; // (expg-1)*expg*( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-2}/mu_{g}

        dlag_dT = pow0*dnug_dT; // derivatives given by the viscosity function
        d2lag_dT2 = pow0*d2nug_dT2; //

        d2lag_dswdT = -expg * pow1*dnug_dT; //
    }

    /* Evaluate sum of la's and its first and second derivatives */
    /* --------------------------------------------------------  */
    la = law + lag;
    dla_dsw = dlaw_dsw + dlag_dsw;
    d2la_dsw2 = d2law_dsw2 + d2lag_dsw2;
    d2la_dswdT = d2law_dswdT + d2lag_dswdT;
    dla_dT = dlaw_dT + dlag_dT;
    d2la_dT2 = d2law_dT2 + d2lag_dT2;

    d2la_dswdT = d2law_dswdT + d2lag_dswdT;

    // Evaluate the Flux Function
    // tw = (vel + lag*(grw - grg));

    flw = law / la;

    // F
    if (degree >= 0) {

        m.set(0, (1. - flw)); // We are giving the value of flg=1-flw, the fractional flow for gas, as required in Flux2Comp2PhasesAdimensionalized as
        // the output of the jet!!! IMPORTANT!

        // Jacobian
        if (degree >= 1) {
            /* Evaluate first derivatives of fluxes relative to gas saturation */
            /* ------------------------------------------------------- */
            dflw_dsw = (la * dlaw_dsw - law * dla_dsw) / (la * la);
            m.set(0, 0, dflw_dsw); /* Now observe that from the chain rule dflg_dsg=dflw_dsw HERE WE DO NOT NEED TO WORRY ABOUT THE DIMENSION */

            /* Evaluate first derivatives of fluxes relative to temperature */
            dflw_dT = (la * dlaw_dT - law * dla_dT) / (la * la);
            m.set(0, 1, (-dflw_dT) * fluxComplete_->TD->T_typical()); // We are giving dflg_dTheta = (  -dflw_dT( T (Theta) )  )*T_typical_
            // Hessian
            if (degree == 2) {
                /* Evaluate second derivatives of fluxes */
                /* ------------------------------------- */
                d2flw_dsw2 = ((la * d2law_dsw2 - law * d2la_dsw2) / la - 2. * dla_dsw * dflw_dsw) / la;

                d2flw_dT2 = (la * la * (la * d2law_dT2 - law * d2la_dT2) - 2. * la * dla_dT * (la * dlaw_dT - law * dla_dT)) / (la * la * la * la);

                /* Mixed derivative */
                /*
                d2flw_dswdT=d2flg_dswdT ;
                 */
                d2flw_dswdT = (la * la * (dla_dT * dlaw_dsw + la * d2law_dswdT - dlaw_dT * dla_dsw - law * d2la_dswdT)
                        - 2. * la * dla_dT * (la * dlaw_dsw - law * dla_dsw)) / (la * la * la * la);

                //               m(0, 0, 0, j000);
                //               m(0, 0, 1, j001);
                //               m(0, 1, 0, j010);
                //               m(0, 1, 1, j011);

                m.set(0, 0, 0, -d2flw_dsw2); // We are passing d2fg_dsg2  NO CHANGE!
                m.set(0, 0, 1, (d2flw_dswdT) * fluxComplete_->TD->T_typical()); // double j001 = d2flg_dsgdTheta = d2flg_dsgdT(T(Theta))*T_typical_  = d2flg_dswdT*T_typical_ = d2flw_dswdT*T_typical_;
                m.set(0, 1, 0, (d2flw_dswdT) * fluxComplete_->TD->T_typical()); // double j010 = j001;
                m.set(0, 1, 1, (-d2flw_dT2) * fluxComplete_->TD->T_typical() * fluxComplete_->TD->T_typical()); // double j011 = d2flg_dTheta2 = d2flg_dT2*T_typical_^2 = -d2flw_dsw2*T_typical_^2 // CHECKED BY HELMUT 21/10/2010

            }
        }
    }
    return 2;
}

int Flux2Comp2PhasesAdimensionalized::FracFlow2PhasesVerticalAdimensionalized::Diff_FracFlow2PhasesVerticalAdimensionalized(double sw, double Theta, int degree, JetMatrix &m) {
    double T = fluxComplete_->TD->Theta2T(Theta);

    double d2Z_dsw2; // <- Unique to Vertical (only different variable with respect to Horizontal)
    double sg, law, lag;
    double dlaw_dsw, dlag_dsw;
    double d2law_dsw2, d2lag_dsw2;
    double la, dla_dsw, d2la_dsw2;
    double dflw_dsw, dflw_dT; // Derivatives of the Fractional Flow of water
    double d2flw_dsw2; // Second derivative of Fractional Flow of water

    double flw; // Fractional Flow for water

    // Some variables to facilitate the operations with powers
    double pow0, pow1, pow2;

    double nuw, nug; // nuw= 1/muw , nug=1/mug
    double dnuw_dT, d2nuw_dT2; // This values are assigned using TD
    // dnuw_dT = dnuw_dT(T).
    // d2nuw_dT2 = d2nuw_dT2(T).

    double dnug_dT, d2nug_dT2; // This values are assigned using TD
    // dnug_dT = dnug_dT(T),
    // d2nug_dT2 = d2nug_dT2(T).

    double dlaw_dT, d2law_dT2; // Temperature derivatives of law (are not passed as parameters).

    double dlag_dT, d2lag_dT2; // Temperature derivatives of lag (are not passed as parameters).

    double d2law_dswdT; // Mixed derivative (it is not passed as parameter).
    double d2lag_dswdT; // Mixed derivative (it is not passed as parameter).

    double dla_dT;
    double d2flw_dswdT;
    double d2la_dswdT;
    double d2flw_dT2;
    double d2la_dT2;

    // Finally we assign sg
    sg = 1. - sw;

    fluxComplete_->TD->inv_muw(T, nuw, dnuw_dT, d2nuw_dT2);
    fluxComplete_->TD->inv_mug(T, nug, dnug_dT, d2nug_dT2);

    double cnw  = fluxComplete_->cnw_parameter->value();
    double cng  = fluxComplete_->cng_parameter->value();
    double expw = fluxComplete_->expw_parameter->value();
    double expg = fluxComplete_->expg_parameter->value();

    // Before evaluating proper
    if (sw < cnw) {
        law = 0.;
        dlaw_dsw = 0.;
        d2law_dsw2 = 0.;

        dlaw_dT = 0.;
        d2law_dT2 = 0.;

        d2law_dswdT = 0.; // remember that d2law_dswdT is equal to d2law_dTdsw
    } else {
        double temp = (sw - cnw) / (1. - cnw - cng);

        pow2 = 0.5 * pow(temp, expw - 2.); // 0.5*( (sw-s_{connate water})/(1-s_{connate water}-s_{con gas}) )^{expw-2}
        // This model coincides with the formulae in Helmut's Thesis. 0.5 is the re-scaling factor.

        pow1 = temp*pow2; // 0.5*( (sw-s_{connate water})/(1-s_{connate water}-s_{con gas}) )^{expw-1}
        // Observe that pow1 = pow((sw - cnw)/(1.-cnw-cng),expw-1) ;

        pow0 = temp*pow1; // 0.5*( (sw-s_{connate water})/(1-s_{connate water-s_{con gas}}) )^{expw}
        // Observe that pow0 = pow((sw - cnw)/(1.-cnw-cng),expw)


        law = pow0*nuw; // k_{rw}/mu_{w}(T) This is equivalent to k_{ra}/mu_{a} in Helmut's thesis.
        dlaw_dsw = expw * pow1*nuw; // expw*( (sw-s_{connate water })/(1.-cnw-cng) )^{expw-1}/mu_{w}  derivation w.r.t sw
        d2law_dsw2 = (expw - 1.) * expw * pow2*nuw; // (expw - 1.) * expw*( (sw-s_{connate water })/(1.-cnw-cng) )^{expw-2}/mu_{w}

        dlaw_dT = pow0*dnuw_dT; // the value of the derivatives of muw
        // will be passed by the viscosity function.
        d2law_dT2 = pow0*d2nuw_dT2; //

        d2law_dswdT = expw * pow1*dnuw_dT;
    }

    if (sg < cng) {
        lag = 0.;
        dlag_dsw = 0.;
        d2lag_dsw2 = 0.;

        dlag_dT = 0.;
        d2lag_dT2 = 0.;

        d2lag_dswdT = 0.;
    } else {

        double temp = (sg - cng) / (1. - cnw - cng);
        pow2 = 0.95 * pow(temp, expg - 2.); // ( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-2}
        pow1 = temp*pow2; // ( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-1} Observe that pow1 = pow(temp, expg-1)
        pow0 = temp*pow1; // ( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg}   Observe that pow0 = pow(temp, expg) ;

        lag = pow0*nug; //  k_{rg}/mu_{g}(T)  this viscosity is function of temperature   ;
        dlag_dsw = -expg * pow1*nug; //  - expg*( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-1}/mu_{g}
        d2lag_dsw2 = (expg - 1.) * expg * pow2*nug; // (expg-1)*expg*( (sg-s_{connate gas})/(1.-cnw-cng) )^{expg-2}/mu_{g}

        dlag_dT = pow0*dnug_dT; // derivatives given by the viscosity function
        d2lag_dT2 = pow0*d2nug_dT2; //

        d2lag_dswdT = -expg * pow1*dnug_dT; //
    }

    /* Evaluate sum of la's and its first and second derivatives */
    /* --------------------------------------------------------  */
    la = law + lag;
    dla_dsw = dlaw_dsw + dlag_dsw;
    d2la_dsw2 = d2law_dsw2 + d2lag_dsw2;
    d2la_dswdT = d2law_dswdT + d2lag_dswdT;
    dla_dT = dlaw_dT + dlag_dT;
    d2la_dT2 = d2law_dT2 + d2lag_dT2;

    d2la_dswdT = d2law_dswdT + d2lag_dswdT;

    flw = law / la;

    // Z
    if (degree >= 0) {

        m.set(0, flw * lag); // double j0 = Z = flw*lag;

        // Jacobian
        if (degree >= 1) {
            /* Evaluate first derivative of Z relative to gas saturation */
            /* --------------------------------------------------------- */

            dflw_dsw = (la * dlaw_dsw - law * dla_dsw) / (la * la); // BE CAREFUL, WE ARE USING THE CHAIN RULE TO PASS THE CORRECT DERIVATIVE
            m.set(0, 0, -(dflw_dsw * lag + flw * dlag_dsw)); // this is double j00 = dZ_dsg= -dZ_dsw = -(dflw_dsw*lag + flw*dlag_dsw); NO CHANGE FOR DIM-LESS

            /* Evaluate first derivative of Z relative to temperature */
            dflw_dT = (la * dlaw_dT - law * dla_dT) / (la * la);

            m.set(0, 1, (dflw_dT * lag + flw * dlag_dT) * fluxComplete_->getThermo()->T_typical()); // this is double j01 =  dZ_dTheta = dZ_dT*T_typical_  = (dflw_dT*lag + flw*dlag_dT)*T_typical_;

            // Hessian
            if (degree == 2) {

                /* Evaluate second derivatives of Z */
                /* ------------------------------------- */
                d2flw_dsw2 = ((la * d2law_dsw2 - law * d2la_dsw2) / la - 2. * dla_dsw * dflw_dsw) / la;

                d2Z_dsw2 = d2flw_dsw2 * lag + 2. * dflw_dsw * dlag_dsw + flw*d2lag_dsw2;

                double j000 = d2Z_dsw2; // Note we are passing d2Z_dsg2 = d2Z_dsw2 NO CHANGE FOR DIM-LESS



                /* Mixed derivative */
                /*
                d2Z_dsgdT = -d2Z_dswdT = -( dflw_dswdT*lag + dflw_dsw*dlag_dT +
                                            dflw_dT*dlag_dsw + flw*d2lag_dswdT  );
                 */


                d2flw_dT2 = (la * la * (la * d2law_dT2 - law * d2la_dT2) - 2. * la * dla_dT * (la * dlaw_dT - law * dla_dT)) / (la * la * la * la);

                d2flw_dswdT = (la * la * (dla_dT * dlaw_dsw + la * d2law_dswdT - dlaw_dT * dla_dsw - law * d2la_dswdT)
                        - 2. * la * dla_dT * (la * dlaw_dsw - law * dla_dsw)) / (la * la * la * la);

                double j001 = (-(d2flw_dswdT * lag + dflw_dsw * dlag_dT +
                        dflw_dT * dlag_dsw + flw * d2lag_dswdT)) * fluxComplete_->getThermo()->T_typical(); // Observe that we are passing really d2Z_dsgdTheta = d2Z_dsgdT*T_typical_

                double j010 = j001;

                /*
                d2Z_dT2 = d2flw_dT2*lag + dflw_dT*dlag_dT +
                          dflw_dT*dlag_dT + flw*d2lag_dT2;
                 */

                double j011 = (d2flw_dT2 * lag + dflw_dT * dlag_dT +
                        dflw_dT * dlag_dT + flw * d2lag_dT2) * fluxComplete_->getThermo()->T_typical() * fluxComplete_->getThermo()->T_typical(); // d2Z_dTheta2 = d2Z_dT2 * T_typical_^2

                m.set(0, 0, 0, j000);
                m.set(0, 0, 1, j001);
                m.set(0, 1, 0, j010);
                m.set(0, 1, 1, j011);
            }
        }
    }
    return 2;
}

void Flux2Comp2PhasesAdimensionalized::type(int t) {
    if (t == FLUX2COMP2PHASESADIMENSIONALIZED_PURE_GRAVITY) {
        has_gravity = true;
        has_horizontal = false;
    } else if (t == FLUX2COMP2PHASESADIMENSIONALIZED_GRAVITY) {
        has_gravity = true;
        has_horizontal = true;
    } else if (t == FLUX2COMP2PHASESADIMENSIONALIZED_HORIZONTAL) {
        has_gravity = false;
        has_horizontal = true;
    }

    return;
}


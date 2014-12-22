#include "Accum2Comp2PhasesAdimensionalized.h"

Accum2Comp2PhasesAdimensionalized::Accum2Comp2PhasesAdimensionalized(Parameter *phi, Thermodynamics *td): AccumulationFunction()
                                                                                                         
{
    phi_parameter_ = phi;
    TD = td;

    reducedAccum_ = new ReducedAccum2Comp2PhasesAdimensionalized(this);
}

Accum2Comp2PhasesAdimensionalized::~Accum2Comp2PhasesAdimensionalized() {
    delete reducedAccum_;
}

// Existe uma discrepancia entre o o significado de s quando este codigo foi
// escrito usando as formulas do Rodrigo, em que s eh a saturacao da fase
// de CO2 supersaturado, e como este codigo vai ser usado no RPN, onde s eh
// a saturacao da agua.
// Isto implica em pequenas mudancas, analogas as que Helmut e eu sugerimos
// para o FracFlow, num email de ontem.
//
// Dan

int Accum2Comp2PhasesAdimensionalized::jet(const WaveState &w, JetMatrix &m, int degree) const {
    double s = w(0); // s_{sigma} = sg in FracFlow2PhasesHorizontal & FracFlow2PhasesVertical
    double Theta = w(1);

    double U = w(2);

    m.resize(3);

    // Some auxiliary variables
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

    // Function
    if (degree >= 0) {
        double phi = phi_parameter_->value();

        double Hr     = Hrj.get(0);
        double Ha     = Haj.get(0);
        double Hsi    = Hsij.get(0);
        double rhosic = rhosicj.get(0);
        double rhosiw = rhosiwj.get(0);
        double rhoac  = rhoacj.get(0);
        double rhoaw  = rhoawj.get(0);

        double out0 = phi * (rhosic * s + rhoac * (1.0 - s)); // G1
        double out1 = phi * (rhosiw * s + rhoaw * (1.0 - s)); // G2
        double out2 = Hr + phi * (Hsi * s + Ha * (1.0 - s)); // G3

//        std::cout << "rhosic = " << rhosic << std::endl;

        m.set(0, out0);
        m.set(1, out1);
        m.set(2, out2);

        // Jacobian
        if (degree >= 1) {
            double d_Hr     = Hrj.get(0, 0);
            double d_Ha     = Haj.get(0, 0);
            double d_Hsi    = Hsij.get(0, 0);
            double d_rhosic = rhosicj.get(0, 0);
            double d_rhosiw = rhosiwj.get(0, 0);
            double d_rhoac  = rhoacj.get(0, 0);
            double d_rhoaw  = rhoawj.get(0, 0);

            double out00 = phi * (rhosic - rhoac); // dG1_ds
            double out01 = phi * (d_rhosic * s + d_rhoac * (1.0 - s)); // dG1_dTheta
            double out02 = 0.; // dG1_dU

            //            printf("out00 = %g\n", out00);
            //            printf("out01 = %g\n", out01);
            //            printf("out02 = %g\n", out02);

            double out10 = phi * (rhosiw - rhoaw); // dG2_ds
            double out11 = phi * (d_rhosiw * s + d_rhoaw * (1.0 - s)); // dG2_dTheta
            double out12 = 0.; // dG1_dU

            double out20 = phi * (Hsi - Ha); // dG3_ds
            double out21 = d_Hr + phi * (d_Hsi * s + d_Ha * (1.0 - s)); // dG3_dTheta
            double out22 = 0.; // dG3_dU

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
            //                    printf("aa(%d, %d) = %g\n", i, j, m(i, j));
            //                }
            //            }

            // Hessian
            if (degree == 2) {
                double d2_Hr     = Hrj.get(0, 0, 0);
                double d2_Ha     = Haj.get(0, 0, 0);
                double d2_Hsi    = Hsij.get(0, 0, 0);
                double d2_rhosic = rhosicj.get(0, 0, 0);
                double d2_rhosiw = rhosiwj.get(0, 0, 0);
                double d2_rhoac  = rhoacj.get(0, 0, 0);
                double d2_rhoaw  = rhoawj.get(0, 0, 0);

                double out000 = 0.; // d2G1_ds2
                double out001 = phi * (d_rhosic - d_rhoac); // d2G1_dsdTheta
                double out002 = 0.; // d2G1_dsdU
                double out010 = out001; // d2G1_dThetads
                double out011 = phi * (d2_rhosic * s + d2_rhoac * (1.0 - s)); // d2G1_dTheta2
                double out012 = 0.; // d2G1_dThetadU
                double out020 = 0.; // d2G1_dUds
                double out021 = 0.; // d2G1_dUdTheta
                double out022 = 0.; // d2G1_dU2

                double out100 = 0.; // d2G2_ds2
                double out101 = phi * (d_rhosiw - d_rhoaw); // d2G2_dsdTheta
                double out102 = 0.; // d2G2_dsdU
                double out110 = phi * (d_rhosiw - d_rhoaw); // d2G2_dThetads
                double out111 = phi * (d2_rhosiw * s + d2_rhoaw * (1.0 - s)); // d2G2_dTheta2
                double out112 = 0.; // d2G2_dThetadU
                double out120 = 0.; // d2G2_dUds
                double out121 = 0.; // d2G2_dUdTheta
                double out122 = 0.; // d2G2_dU2

                double out200 = 0.; // d2G3_ds2
                double out201 = phi * (d_Hsi - d_Ha); // d2G3_dsdTheta
                double out202 = 0.; // d2G3_dsdU
                double out210 = out201; // d2G3_dThetads
                double out211 = d2_Hr + phi * (d2_Hsi * s + d2_Ha * (1.0 - s)); // d2G3_dTheta2
                double out212 = 0.; // d2G3_dThetadU
                double out220 = 0.; // d2G3_dUds
                double out221 = 0.; // d2G3_dUdTheta
                double out222 = 0.; // d2G3_dU2

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


Accum2Comp2PhasesAdimensionalized::ReducedAccum2Comp2PhasesAdimensionalized::ReducedAccum2Comp2PhasesAdimensionalized(Accum2Comp2PhasesAdimensionalized * outer):AccumulationFunction(outer->accumulationParams()){
    totalAccum_=outer;

 }

Accum2Comp2PhasesAdimensionalized::ReducedAccum2Comp2PhasesAdimensionalized::~ReducedAccum2Comp2PhasesAdimensionalized(){

}

 Accum2Comp2PhasesAdimensionalized::ReducedAccum2Comp2PhasesAdimensionalized * Accum2Comp2PhasesAdimensionalized::getReducedAccumulation()const{
     return reducedAccum_;

 }

int Accum2Comp2PhasesAdimensionalized::ReducedAccum2Comp2PhasesAdimensionalized::jet(const WaveState &w, JetMatrix &m, int degree) const{
    double s     = w(0); // s_{sigma} = sg in FracFlow2PhasesHorizontal & FracFlow2PhasesVertical
    double Theta = w(1);

    m.resize(2, 3);

    // Some auxiliary variables
    JetMatrix Hrj(1);
    totalAccum_->TD->RockEnthalpyVol_jet(Theta, degree, Hrj);

    JetMatrix Haj(1);
    totalAccum_->TD->AqueousEnthalpyVol_jet(Theta, degree, Haj);

    JetMatrix Hsij(1);
    totalAccum_->TD->SuperCriticEnthalpyVol_jet(Theta, degree, Hsij);

    JetMatrix rhosicj(1);
    totalAccum_->TD->Rhosic_jet(Theta, degree, rhosicj);

    JetMatrix rhosiwj(1);
    totalAccum_->TD->Rhosiw_jet(Theta, degree, rhosiwj);

    JetMatrix rhoacj(1);
    totalAccum_->TD->Rhoac_jet(Theta, degree, rhoacj);

    JetMatrix rhoawj(1);
    totalAccum_->TD->Rhoaw_jet(Theta, degree, rhoawj);

    // Function
    if (degree >= 0){
        double phi = totalAccum_->phi_parameter_->value();        

        double Hr     = Hrj.get(0);
        double Ha     = Haj.get(0);
        double Hsi    = Hsij.get(0);
        double rhosic = rhosicj.get(0);
        double rhosiw = rhosiwj.get(0);
        double rhoac  = rhoacj.get(0);
        double rhoaw  = rhoawj.get(0);

        double out0 = phi*(rhosic*s + rhoac*(1.0 - s)); // G1
        double out1 = phi*(rhosiw*s + rhoaw*(1.0 - s)); // G2
        double out2 = Hr + phi*(Hsi*s + Ha*(1.0 - s) ); // G3

        m.set(0, out0);
        m.set(1, out1);
        m.set(2, out2);

        // Jacobian
        if (degree >= 1){
            double d_Hr     = Hrj.get(0, 0);
            double d_Ha     = Haj.get(0, 0);
            double d_Hsi    = Hsij.get(0, 0);
            double d_rhosic = rhosicj.get(0, 0);
            double d_rhosiw = rhosiwj.get(0, 0);
            double d_rhoac  = rhoacj.get(0, 0);
            double d_rhoaw  = rhoawj.get(0, 0);

            double out00 = phi*(rhosic - rhoac); // dG1_ds
            double out01 = phi*(d_rhosic*s + d_rhoac*(1.0 - s)); // dG1_dTheta

//            printf("out00 = %g\n", out00);
//            printf("out01 = %g\n", out01);

            double out10 = phi*(rhosiw - rhoaw); // dG2_ds
            double out11 = phi*(d_rhosiw*s + d_rhoaw*(1.0 - s)); // dG2_dTheta

            double out20 = phi*(Hsi-Ha); // dG3_ds
            double out21 = d_Hr + phi*( d_Hsi*s + d_Ha*(1.0 - s) ); // dG3_dTheta

            m.set(0, 0, out00);
            m.set(0, 1, out01);

            m.set(1, 0, out10);
            m.set(1, 1, out11);

            m.set(2, 0, out20);
            m.set(2, 1, out21);

//            for (int i = 0; i < 3; i++){
//                for (int j = 0; j < 3; j++){
//                    printf("aa(%d, %d) = %g\n", i, j, m(i, j));
//                }
//            }


            // Hessian
            if (degree == 2){
                double d2_Hr     = Hrj.get(0, 0, 0);
                double d2_Ha     = Haj.get(0, 0, 0);
                double d2_Hsi    = Hsij.get(0, 0, 0);
                double d2_rhosic = rhosicj.get(0, 0, 0);
                double d2_rhosiw = rhosiwj.get(0, 0, 0);
                double d2_rhoac  = rhoacj.get(0, 0, 0);
                double d2_rhoaw  = rhoawj.get(0, 0, 0);

                double out000 = 0.; // d2G1_ds2
                double out001 = phi*(d_rhosic - d_rhoac); // d2G1_dsdTheta
                double out010 = out001; // d2G1_dThetads
                double out011 = phi*(d2_rhosic*s + d2_rhoac*(1.0 - s)); // d2G1_dTheta2

                double out100 = 0.; // d2G2_ds2
                double out101 = phi*(d_rhosiw - d_rhoaw); // d2G2_dsdTheta
                double out110 = phi*(d_rhosiw - d_rhoaw); // d2G2_dThetads
                double out111 = phi*(d2_rhosiw*s + d2_rhoaw*(1.0 - s)); // d2G2_dTheta2

                double out200 = 0.; // d2G3_ds2
                double out201 = phi*(d_Hsi-d_Ha); // d2G3_dsdTheta
                double out210 = out201; // d2G3_dThetads
                double out211 = d2_Hr + phi*(d2_Hsi*s + d2_Ha*(1.0 - s) ); // d2G3_dTheta2

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

    return  2; //SUCCESSFUL_PROCEDURE;
}

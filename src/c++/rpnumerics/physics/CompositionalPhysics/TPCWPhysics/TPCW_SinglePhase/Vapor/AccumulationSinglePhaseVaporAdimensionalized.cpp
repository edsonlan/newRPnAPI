#include "AccumulationSinglePhaseVaporAdimensionalized.h"

//AccumulationSinglePhaseVaporAdimensionalized::AccumulationSinglePhaseVaporAdimensionalized(const AccumulationSinglePhaseVaporAdimensionalized &copy) : AccumulationFunction(copy.accumulationParams()) {

//    thermo = copy.thermo;
//}


AccumulationSinglePhaseVaporAdimensionalized::AccumulationSinglePhaseVaporAdimensionalized(double phi, Thermodynamics *t):AccumulationFunction(AccumulationParams()) {
//    accumulationParams().component(0,phi);

    phi_ = phi;
    thermo = t;
}

AccumulationSinglePhaseVaporAdimensionalized * AccumulationSinglePhaseVaporAdimensionalized::clone() const {
    return new AccumulationSinglePhaseVaporAdimensionalized(*this);
}

AccumulationSinglePhaseVaporAdimensionalized::~AccumulationSinglePhaseVaporAdimensionalized() {
    //    delete thermo;
}

int AccumulationSinglePhaseVaporAdimensionalized::jet(const WaveState &w, JetMatrix &m, int degree) const {
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE



//    double phi_ = accumulationParams().component(0);


    double yw = w(0); // s_{sigma} = sg in FracFlow2PhasesHorizontal & FracFlow2PhasesVertical
    double Theta = w(1);
    double U = w(2);

    // Some auxiliary variables
    JetMatrix Hsij(2);
    thermo->SuperCriticEnthalpyVol_jet(yw, Theta, degree, Hsij);

    JetMatrix rhosicj(2);
    thermo->Rhosic_jet(yw, Theta, degree, rhosicj);

    JetMatrix rhosiwj(2);
    thermo->Rhosiw_jet(yw, Theta, degree, rhosiwj);

    JetMatrix revj(1);
    thermo->RockEnthalpyVol_jet(Theta, degree, revj);

    if (degree >= 0) {
        double Hsi = Hsij.get(0);
        double rhosic = rhosicj.get(0);
        double rhosiw = rhosiwj.get(0);
        double Hr = revj.get(0);

        double G0 = phi_*rhosic;
        double G1 = phi_*rhosiw;
        double G2 = Hr + phi_*Hsi;

        m.set(0, G0);
        m.set(1, G1);
        m.set(2, G2);

        if (degree >= 1) {
            double dHsi_dyw = Hsij.get(0, 0);
            double dHsi_dTheta = Hsij.get(0, 1);

            double drhosic_dyw = rhosicj.get(0, 0);
            double drhosic_dTheta = rhosicj.get(0, 1);

            double drhosiw_dyw = rhosiwj.get(0, 0);
            double drhosiw_dTheta = rhosiwj.get(0, 1);

            double dHr_dTheta = revj.get(0, 0);

            double dG0_dyw = phi_*drhosic_dyw;
            double dG0_dTheta = phi_*drhosic_dTheta;
            double dG0_dU = 0.0;

            double dG1_dyw = phi_*drhosiw_dyw;
            double dG1_dTheta = phi_*drhosiw_dTheta;
            double dG1_dU = 0.0;

            double dG2_dyw = phi_*dHsi_dyw;
            double dG2_dTheta = dHr_dTheta + phi_*dHsi_dTheta;
            double dG2_dU = 0.0;

            m.set(0, 0, dG0_dyw);
            m.set(0, 1, dG0_dTheta);
            m.set(0, 2, dG0_dU);

            m.set(1, 0, dG1_dyw);
            m.set(1, 1, dG1_dTheta);
            m.set(1, 2, dG1_dU);

            m.set(2, 0, dG2_dyw);
            m.set(2, 1, dG2_dTheta);
            m.set(2, 2, dG2_dU);

            if (degree >= 2) {
                double d2Hsi_dyw2 = Hsij.get(0, 0, 0);
                double d2Hsi_dywdTheta = Hsij.get(0, 0, 1);
                double d2Hsi_dThetadyw = Hsij.get(0, 1, 0);
                double d2Hsi_dTheta2 = Hsij.get(0, 1, 1);

                double d2rhosic_dyw2 = rhosicj.get(0, 0, 0);
                double d2rhosic_dywdTheta = rhosicj.get(0, 0, 1);
                double d2rhosic_dThetadyw = rhosicj.get(0, 1, 0);
                double d2rhosic_dTheta2 = rhosicj.get(0, 1, 1);

                double d2rhosiw_dyw2 = rhosiwj.get(0, 0, 0);
                double d2rhosiw_dywdTheta = rhosiwj.get(0, 0, 1);
                double d2rhosiw_dThetadyw = rhosiwj.get(0, 1, 0);
                double d2rhosiw_dTheta2 = rhosiwj.get(0, 1, 1);

                double d2Hr_dTheta2 = revj.get(0, 0, 0);

                double d2G0dyw2 = phi_*d2rhosic_dyw2;
                double d2G0dywdTheta = phi_*d2rhosic_dywdTheta;
                double d2G0dywdU = 0.0;

                double d2G0dThetadyw = d2G0dywdTheta;
                double d2G0dTheta2 = phi_*d2rhosic_dTheta2;
                double d2G0dThetadU = 0.0;

                double d2G0dUdyw = d2G0dywdU;
                double d2G0dUdTheta = d2G0dThetadU;
                double d2G0dU2 = 0.0;

                m.set(0, 0, 0, d2G0dyw2);
                m.set(0, 0, 1, d2G0dywdTheta);
                m.set(0, 0, 2, d2G0dywdU);
                m.set(0, 1, 0, d2G0dThetadyw);
                m.set(0, 1, 1, d2G0dTheta2);
                m.set(0, 1, 2, d2G0dThetadU);
                m.set(0, 2, 0, d2G0dUdyw);
                m.set(0, 2, 1, d2G0dUdTheta);
                m.set(0, 2, 2, d2G0dU2);

                double d2G1dyw2 = phi_*d2rhosiw_dyw2;
                double d2G1dywdTheta = phi_*d2rhosiw_dywdTheta;
                double d2G1dywdU = 0.0;

                double d2G1dThetadyw = d2G1dywdTheta;
                double d2G1dTheta2 = phi_*d2rhosiw_dTheta2;
                double d2G1dThetadU = 0.0;

                double d2G1dUdyw = d2G1dywdU;
                double d2G1dUdTheta = d2G1dThetadU;
                double d2G1dU2 = 0.0;

                m.set(1, 0, 0, d2G1dyw2);
                m.set(1, 0, 1, d2G1dywdTheta);
                m.set(1, 0, 2, d2G1dywdU);
                m.set(1, 1, 0, d2G1dThetadyw);
                m.set(1, 1, 1, d2G1dTheta2);
                m.set(1, 1, 2, d2G1dThetadU);
                m.set(1, 2, 0, d2G1dUdyw);
                m.set(1, 2, 1, d2G1dUdTheta);
                m.set(1, 2, 2, d2G1dU2);

                double d2G2dyw2 = phi_*d2Hsi_dyw2;
                double d2G2dywdTheta = phi_*d2Hsi_dywdTheta;
                double d2G2dywdU = 0.0;

                double d2G2dThetadyw = d2G2dywdTheta;
                double d2G2dTheta2 = phi_ * d2Hsi_dTheta2 + d2Hr_dTheta2;
                double d2G2dThetadU = 0.0;

                double d2G2dUdyw = d2G2dywdU;
                double d2G2dUdTheta = d2G2dThetadU;
                double d2G2dU2 = 0.0;

                m.set(2, 0, 0, d2G2dyw2);
                m.set(2, 0, 1, d2G2dywdTheta);
                m.set(2, 0, 2, d2G2dywdU);
                m.set(2, 1, 0, d2G2dThetadyw);
                m.set(2, 1, 1, d2G2dTheta2);
                m.set(2, 1, 2, d2G2dThetadU);
                m.set(2, 2, 0, d2G2dUdyw);
                m.set(2, 2, 1, d2G2dUdTheta);
                m.set(2, 2, 2, d2G2dU2);
            }
        }
    }

    return 2; //SUCCESSFUL_PROCEDURE;
}


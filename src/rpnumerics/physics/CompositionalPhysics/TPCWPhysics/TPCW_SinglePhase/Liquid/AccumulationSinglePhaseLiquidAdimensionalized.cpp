#include "AccumulationSinglePhaseLiquidAdimensionalized.h"

//AccumulationSinglePhaseLiquidAdimensionalized::AccumulationSinglePhaseLiquidAdimensionalized(const AccumulationSinglePhaseLiquidAdimensionalized &copy) : AccumulationFunction(copy.accumulationParams()){

//    thermo = copy.thermo;
//}

AccumulationSinglePhaseLiquidAdimensionalized::AccumulationSinglePhaseLiquidAdimensionalized(double phi, Thermodynamics *t) {

//    accumulationParams().component(0,phi);
    phi_ = phi;
    thermo = t;
}

AccumulationSinglePhaseLiquidAdimensionalized * AccumulationSinglePhaseLiquidAdimensionalized::clone() const {
    return new AccumulationSinglePhaseLiquidAdimensionalized(*this);
}

AccumulationSinglePhaseLiquidAdimensionalized::~AccumulationSinglePhaseLiquidAdimensionalized(){
//    delete thermo;
}

int AccumulationSinglePhaseLiquidAdimensionalized::jet(const WaveState &w, JetMatrix &m, int degree) const{
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE
    
    
//    double phi_=accumulationParams().component(0);

    double xc    = w(0); // s_{sigma} = sg in FracFlow2PhasesHorizontal & FracFlow2PhasesVertical
    double Theta = w(1);
    double U     = w(2);

    // Some auxiliary variables
    JetMatrix Haj(2);
    thermo->AqueousEnthalpyVol_jet(xc, Theta, degree, Haj);

    JetMatrix rhoacj(2);
    thermo->Rhoac_jet(xc, Theta, degree, rhoacj);

    JetMatrix rhoawj(2);
    thermo->Rhoaw_jet(xc, Theta, degree, rhoawj);

    JetMatrix revj(1);
    thermo->RockEnthalpyVol_jet(Theta, degree, revj);

    if (degree >= 0){
        double Ha    = Haj.get(0);
        double rhoac = rhoacj.get(0);
        double rhoaw = rhoawj.get(0);
        double Hr    = revj.get(0);

        double G0 = phi_*rhoac;
        double G1 = phi_*rhoaw;
        double G2 = Hr + phi_*Ha;

        m.set(0, G0);
        m.set(1, G1);
        m.set(2, G2);

        if (degree >= 1){
            double dHa_dxc       = Haj.get(0, 0);
            double dHa_dTheta    = Haj.get(0, 1);

            double drhoac_dxc    = rhoacj.get(0, 0);
            double drhoac_dTheta = rhoacj.get(0, 1);

            double drhoaw_dxc    = rhoawj.get(0, 0);
            double drhoaw_dTheta = rhoawj.get(0, 1);

            double dHr_dTheta     = revj.get(0, 0);

            double dG0_dxc        = phi_*drhoac_dxc;
            double dG0_dTheta     = phi_*drhoac_dTheta;
            double dG0_dU         = 0.0;

            double dG1_dxc        = phi_*drhoaw_dxc;
            double dG1_dTheta     = phi_*drhoaw_dTheta;
            double dG1_dU         = 0.0;

            double dG2_dxc        = phi_*dHa_dxc;
            double dG2_dTheta     = dHr_dTheta + phi_*dHa_dTheta;
            double dG2_dU         = 0.0;

            m.set(0, 0, dG0_dxc);
            m.set(0, 1, dG0_dTheta);
            m.set(0, 2, dG0_dU);

            m.set(1, 0, dG1_dxc);
            m.set(1, 1, dG1_dTheta);
            m.set(1, 2, dG1_dU);

            m.set(2, 0, dG2_dxc);
            m.set(2, 1, dG2_dTheta);
            m.set(2, 2, dG2_dU);

            if (degree >= 2){
                double d2Ha_dxc2         = Haj.get(0, 0, 0);
                double d2Ha_dxcdTheta    = Haj.get(0, 0, 1);
                double d2Ha_dThetadxc    = Haj.get(0, 1, 0);
                double d2Ha_dTheta2      = Haj.get(0, 1, 1);

                double d2rhoac_dxc2      = rhoacj.get(0, 0, 0);
                double d2rhoac_dxcdTheta = rhoacj.get(0, 0, 1);
                double d2rhoac_dThetadxc = rhoacj.get(0, 1, 0);
                double d2rhoac_dTheta2   = rhoacj.get(0, 1, 1);

                double d2rhoaw_dxc2      = rhoawj.get(0, 0, 0);
                double d2rhoaw_dxcdTheta = rhoawj.get(0, 0, 1);
                double d2rhoaw_dThetadxc = rhoawj.get(0, 1, 0);
                double d2rhoaw_dTheta2   = rhoawj.get(0, 1, 1);

                double d2Hr_dTheta2       = revj.get(0, 0, 0);

                double d2G0dxc2      = phi_*d2rhoac_dxc2;
                double d2G0dxcdTheta = phi_*d2rhoac_dxcdTheta;
                double d2G0dxcdU     = 0.0;

                double d2G0dThetadxc = d2G0dxcdTheta;
                double d2G0dTheta2   = phi_*d2rhoac_dTheta2;
                double d2G0dThetadU  = 0.0;

                double d2G0dUdxc     = d2G0dxcdU;
                double d2G0dUdTheta  = d2G0dThetadU;
                double d2G0dU2       = 0.0;

                m.set(0, 0, 0, d2G0dxc2);
                m.set(0, 0, 1, d2G0dxcdTheta);
                m.set(0, 0, 2, d2G0dxcdU);
                m.set(0, 1, 0, d2G0dThetadxc);
                m.set(0, 1, 1, d2G0dTheta2);
                m.set(0, 1, 2, d2G0dThetadU);
                m.set(0, 2, 0, d2G0dUdxc);
                m.set(0, 2, 1, d2G0dUdTheta);
                m.set(0, 2, 2, d2G0dU2);

                double d2G1dxc2      = phi_*d2rhoaw_dxc2;
                double d2G1dxcdTheta = phi_*d2rhoaw_dxcdTheta;
                double d2G1dxcdU     = 0.0;

                double d2G1dThetadxc = d2G1dxcdTheta;
                double d2G1dTheta2   = phi_*d2rhoaw_dTheta2;
                double d2G1dThetadU  = 0.0;

                double d2G1dUdxc     = d2G1dxcdU;
                double d2G1dUdTheta  = d2G1dThetadU;
                double d2G1dU2       = 0.0;

                m.set(1, 0, 0, d2G1dxc2);
                m.set(1, 0, 1, d2G1dxcdTheta);
                m.set(1, 0, 2, d2G1dxcdU);
                m.set(1, 1, 0, d2G1dThetadxc);
                m.set(1, 1, 1, d2G1dTheta2);
                m.set(1, 1, 2, d2G1dThetadU);
                m.set(1, 2, 0, d2G1dUdxc);
                m.set(1, 2, 1, d2G1dUdTheta);
                m.set(1, 2, 2, d2G1dU2);

                double d2G2dxc2      = phi_*d2Ha_dxc2;
                double d2G2dxcdTheta = phi_*d2Ha_dxcdTheta;
                double d2G2dxcdU     = 0.0;

                double d2G2dThetadxc = d2G2dxcdTheta;
                double d2G2dTheta2   = phi_*d2Ha_dTheta2 + d2Hr_dTheta2;
                double d2G2dThetadU  = 0.0;

                double d2G2dUdxc     = d2G2dxcdU;
                double d2G2dUdTheta  = d2G2dThetadU;
                double d2G2dU2       = 0.0;

                m.set(2, 0, 0, d2G2dxc2);
                m.set(2, 0, 1, d2G2dxcdTheta);
                m.set(2, 0, 2, d2G2dxcdU);
                m.set(2, 1, 0, d2G2dThetadxc);
                m.set(2, 1, 1, d2G2dTheta2);
                m.set(2, 1, 2, d2G2dThetadU);
                m.set(2, 2, 0, d2G2dUdxc);
                m.set(2, 2, 1, d2G2dUdTheta);
                m.set(2, 2, 2, d2G2dU2);
            }
        }
    }

    return  2; //SUCCESSFUL_PROCEDURE;
}


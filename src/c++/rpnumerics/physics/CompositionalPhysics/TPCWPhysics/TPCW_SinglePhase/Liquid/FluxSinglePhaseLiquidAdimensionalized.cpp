#include "FluxSinglePhaseLiquidAdimensionalized.h"

FluxSinglePhaseLiquidAdimensionalized::FluxSinglePhaseLiquidAdimensionalized(const FluxSinglePhaseLiquidAdimensionalized &copy):FluxFunction(copy.fluxParams()) {
   
    thermo = copy.thermo;
}

FluxSinglePhaseLiquidAdimensionalized::FluxSinglePhaseLiquidAdimensionalized(Thermodynamics *t):FluxFunction(FluxParams(RealVector(2))){

    thermo = t;
}

FluxSinglePhaseLiquidAdimensionalized * FluxSinglePhaseLiquidAdimensionalized::clone() const {
    return new FluxSinglePhaseLiquidAdimensionalized(*this);
}

FluxSinglePhaseLiquidAdimensionalized::~FluxSinglePhaseLiquidAdimensionalized(){
    //delete thermo;
}

int FluxSinglePhaseLiquidAdimensionalized::jet(const WaveState &w, JetMatrix &m, int degree) const{
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE

    double xc    = w(0); // s_{sigma} = sg in FracFlow2PhasesHorizontal & FracFlow2PhasesVertical
    double Theta = w(1);
    double U     = w(2);

    

//    printf("xc = %g; Theta = %g; U = %g\n", xc, Theta, U);

    // Some auxiliary variables
    JetMatrix Haj(2);
    thermo->AqueousEnthalpyVol_jet(xc, Theta, degree, Haj); //printf("Done: AqueousEnthalpyVol_jet\n");

    JetMatrix rhoacj(2);
    thermo->Rhoac_jet(xc, Theta, degree, rhoacj); //printf("Done: Rhoac_jet\n");

    JetMatrix rhoawj(2);
    thermo->Rhoaw_jet(xc, Theta, degree, rhoawj); //printf("Done: Rhoaw_jet\n");

    if (degree >= 0){
        double Ha    = Haj.get(0);
        double rhoac = rhoacj.get(0);
        double rhoaw = rhoawj.get(0);

        double F0 = U*rhoac;
        double F1 = U*rhoaw;
        double F2 = U*Ha;

        m.set(0, F0);
        m.set(1, F1);
        m.set(2, F2);

        if (degree >= 1){
            double dHa_dxc       = Haj.get(0, 0);
            double dHa_dTheta    = Haj.get(0, 1);

            double drhoac_dxc    = rhoacj.get(0, 0);
            double drhoac_dTheta = rhoacj.get(0, 1);

            double drhoaw_dxc    = rhoawj.get(0, 0);
            double drhoaw_dTheta = rhoawj.get(0, 1);

            double dF0_dxc        = U*drhoac_dxc;
            double dF0_dTheta     = U*drhoac_dTheta;
            double dF0_dU         = rhoac;

    
            double dF1_dxc        = U*drhoaw_dxc;
            double dF1_dTheta     = U*drhoaw_dTheta;
            double dF1_dU         = rhoaw;

            double dF2_dxc        = U*dHa_dxc;
            double dF2_dTheta     = U*dHa_dTheta;
            double dF2_dU         = Ha;

            m.set(0, 0, dF0_dxc);
            m.set(0, 1, dF0_dTheta);
            m.set(0, 2, dF0_dU);

            m.set(1, 0, dF1_dxc);
            m.set(1, 1, dF1_dTheta);
            m.set(1, 2, dF1_dU);

            m.set(2, 0, dF2_dxc);
            m.set(2, 1, dF2_dTheta);
            m.set(2, 2, dF2_dU);

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

                double d2F0dxc2      = U*d2rhoac_dxc2;
                double d2F0dxcdTheta = U*d2rhoac_dxcdTheta;
                double d2F0dxcdU     = drhoac_dxc;

                double d2F0dThetadxc = d2F0dxcdTheta;
                double d2F0dTheta2   = U*d2rhoac_dTheta2;
                double d2F0dThetadU  = drhoac_dTheta;

                double d2F0dUdxc     = d2F0dxcdU;
                double d2F0dUdTheta  = d2F0dThetadU;
                double d2F0dU2       = 0.0;

                m.set(0, 0, 0, d2F0dxc2);
                m.set(0, 0, 1, d2F0dxcdTheta);
                m.set(0, 0, 2, d2F0dxcdU);
                m.set(0, 1, 0, d2F0dThetadxc);
                m.set(0, 1, 1, d2F0dTheta2);
                m.set(0, 1, 2, d2F0dThetadU);
                m.set(0, 2, 0, d2F0dUdxc);
                m.set(0, 2, 1, d2F0dUdTheta);
                m.set(0, 2, 2, d2F0dU2);

                double d2F1dxc2      = U*d2rhoaw_dxc2;
                double d2F1dxcdTheta = U*d2rhoaw_dxcdTheta;
                double d2F1dxcdU     = drhoaw_dxc;

                double d2F1dThetadxc = d2F1dxcdTheta;
                double d2F1dTheta2   = U*d2rhoaw_dTheta2;
                double d2F1dThetadU  = drhoaw_dTheta;

                double d2F1dUdxc     = d2F1dxcdU;
                double d2F1dUdTheta  = d2F1dThetadU;
                double d2F1dU2       = 0.0;

                m.set(1, 0, 0, d2F1dxc2);
                m.set(1, 0, 1, d2F1dxcdTheta);
                m.set(1, 0, 2, d2F1dxcdU);
                m.set(1, 1, 0, d2F1dThetadxc);
                m.set(1, 1, 1, d2F1dTheta2);
                m.set(1, 1, 2, d2F1dThetadU);
                m.set(1, 2, 0, d2F1dUdxc);
                m.set(1, 2, 1, d2F1dUdTheta);
                m.set(1, 2, 2, d2F1dU2);

                double d2F2dxc2      = U*d2Ha_dxc2;
                double d2F2dxcdTheta = U*d2Ha_dxcdTheta;
                double d2F2dxcdU     = dHa_dxc;

                double d2F2dThetadxc = d2F2dxcdTheta;
                double d2F2dTheta2   = U*d2Ha_dTheta2;
                double d2F2dThetadU  = dHa_dTheta;

                double d2F2dUdxc     = d2F2dxcdU;
                double d2F2dUdTheta  = d2F2dThetadU;
                double d2F2dU2       = 0.0;

                m.set(2, 0, 0, d2F2dxc2);
                m.set(2, 0, 1, d2F2dxcdTheta);
                m.set(2, 0, 2, d2F2dxcdU);
                m.set(2, 1, 0, d2F2dThetadxc);
                m.set(2, 1, 1, d2F2dTheta2);
                m.set(2, 1, 2, d2F2dThetadU);
                m.set(2, 2, 0, d2F2dUdxc);
                m.set(2, 2, 1, d2F2dUdTheta);
                m.set(2, 2, 2, d2F2dU2);
            }
        }
    }

    return  2; //SUCCESSFUL_PROCEDURE;
}


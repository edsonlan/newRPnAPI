#include "DeadVolatileVolatileGasHydrodynamics.h"

DeadVolatileVolatileGasHydrodynamics::DeadVolatileVolatileGasHydrodynamics(DeadVolatileVolatileGasThermodynamics *t){
    thermo = t;
    
    info_auxiliaryfunction_ = std::string("Hydrodynamics");
}

DeadVolatileVolatileGasHydrodynamics::~DeadVolatileVolatileGasHydrodynamics(){
}

// r is R.get(0).
//
void DeadVolatileVolatileGasHydrodynamics::oil_fractional_flow(int degree, double s, double r, JetMatrix &fo_jet){
    fo_jet.resize(2, 1);

    if (degree >= 0){
        double M = s*s + r*(1.0 - s)*(1.0 - s);
        double inv_M = 1.0/M;

        double fo = s*s*inv_M;
        fo_jet.set(0, fo);

        if (degree >= 1){
            double dM_ds = 2.0*(s - r*(1.0 - s));
            double dM_dr = (1.0 - s)*(1.0 - s);

            double FMs = dM_ds*inv_M;
            double FMr = dM_dr*inv_M;

            double dfo_ds = 2.0*s*inv_M - FMs*fo;
            fo_jet.set(0, 0, dfo_ds);

            double dfo_dr = -FMr*fo;
            fo_jet.set(0, 1, dfo_dr);

            if (degree >= 2){
                double d2M_ds2  = 2.0*(1.0 + r);
                double d2M_dsdr = -2.0*(1.0 - s);
                // double d2M_dr2  = 0.0; // Avoid using this zero.

                double FMss = d2M_ds2*inv_M;
                double FMsr = d2M_dsdr*inv_M;
                double FMrr = 0.0; // Avoid using this zero.

                fo_jet.set(0, 0, 0, 2.0*inv_M - 2.0*FMs*dfo_ds - FMss*fo);
                fo_jet.set(0, 0, 1, -FMs*dfo_dr - FMr*dfo_ds - FMsr*fo);
                fo_jet.set(0, 1, 0, fo_jet.get(0, 0, 1));
                fo_jet.set(0, 1, 1, 2.0*FMr*FMr*fo);
            }
        }
    }

    return;
}

void DeadVolatileVolatileGasHydrodynamics::fractional_flow(int degree, double s, double y, JetMatrix &f_jet){
    f_jet.resize(2, 1);

    JetMatrix R_jet;
    thermo->viscosity_ratio(degree, y, R_jet);

    if (degree >= 0){
        double r = R_jet.get(0);
        JetMatrix fo_jet;
        oil_fractional_flow(degree, s, r, fo_jet);

        f_jet.set(0, fo_jet.get(0)); // f = fo.

        if (degree >= 1){
            f_jet.set(0, 0, fo_jet.get(0, 0)); // df_ds = dfo_ds.

            double dR_dy  = R_jet.get(0, 0);
            double dfo_dr = fo_jet.get(0, 1);
            f_jet.set(0, 1, dfo_dr*dR_dy); // df_dr = dfo_dr*dR_dy.

            if (degree >= 2){
                f_jet.set(0, 0, 0, fo_jet.get(0, 0, 0));                 // d2f_ds2  = d2fo_ds2.
                
                double d2fo_dsdr = fo_jet.get(0, 0, 1);
                f_jet.set(0, 0, 1, d2fo_dsdr*dR_dy); // d2f_dsdr = d2fo_dsdr*dR_dy.
                
                f_jet.set(0, 1, 0, f_jet.get(0, 0, 1));
                
                double d2fo_dr2 = fo_jet.get(0, 1, 1);
                double d2R_dy2 = R_jet.get(0, 0, 0);
                f_jet.set(0, 1, 1, d2fo_dr2*dR_dy*dR_dy + dfo_dr*d2R_dy2); 
            }
        }
    }

    return;
}

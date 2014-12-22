#include "DeadVolatileVolatileGasFluxFunction.h"

DeadVolatileVolatileGasFluxFunction::DeadVolatileVolatileGasFluxFunction(DeadVolatileVolatileGasThermodynamics *th, DeadVolatileVolatileGasHydrodynamics *hy) : FluxFunction() {
    thermo = th;
    hydro  = hy;
}

DeadVolatileVolatileGasFluxFunction::~DeadVolatileVolatileGasFluxFunction(){
}

int DeadVolatileVolatileGasFluxFunction::reduced_jet(const WaveState &state, JetMatrix &m, int degree) const {
    m.resize(2, 3);

    double s = state(0);
    double y = state(1);

    // Compute the molar densities.
    //
    JetMatrix Rga_jet, Rgb_jet;
    thermo->gas_molar_density_a(degree, y, Rga_jet);
    thermo->gas_molar_density_b(degree, y, Rgb_jet);
    
    JetMatrix Ro_jet;
    thermo->oil_molar_density(degree, y, Ro_jet);
    
    JetMatrix Rob_jet, Rod_jet;
    thermo->oil_molar_density_b(degree, Ro_jet, y, Rob_jet);
    thermo->oil_molar_density_d(degree, Ro_jet, y, Rod_jet);
    
    JetMatrix f_jet;
    hydro->fractional_flow(degree, s, y, f_jet);
    
    if (degree >= 0){
        double Rga = Rga_jet.get(0);
        double Rgb = Rgb_jet.get(0);
        double Rob = Rob_jet.get(0);
        double Rod = Rod_jet.get(0);
        
        double f = f_jet.get(0);
        
        m.set(0, Rga*(1.0 - f));
        m.set(1, Rgb + (Rob - Rgb)*f);
        m.set(2, Rod*f);
        
        if (degree >= 1){
            double dRga_dy = Rga_jet.get(0, 0);
            double dRgb_dy = Rgb_jet.get(0, 0);
            
            double dRob_dy = Rob_jet.get(0, 0);
            double dRod_dy = Rod_jet.get(0, 0);
            
            double df_ds = f_jet.get(0, 0);
            double df_dy = f_jet.get(0, 1);
            
            m.set(0, 0, -Rga*df_ds);
            m.set(0, 1, dRga_dy*(1.0 - f) - Rga*df_dy);

            m.set(1, 0, (Rob - Rgb)*df_ds);
            m.set(1, 1, dRgb_dy + (dRob_dy - dRgb_dy)*f + (Rob - Rgb)*df_dy);

            m.set(2, 0, Rod*df_ds);
            m.set(2, 1, dRod_dy*f + Rod*df_dy);
            
            if (degree >= 2){
                double d2Rga_dy2 = Rga_jet.get(0, 0, 0);
                double d2Rgb_dy2 = Rgb_jet.get(0, 0, 0);
                double d2Rob_dy2 = Rob_jet.get(0, 0, 0);
                double d2Rod_dy2 = Rod_jet.get(0, 0, 0);

                double d2f_ds2  = f_jet.get(0, 0, 0);
                double d2f_dsdy = f_jet.get(0, 0, 1);
                double d2f_dyds = d2f_dsdy;
                double d2f_dy2  = f_jet.get(0, 1, 1);

                m.set(0, 0, 0, -Rga*d2f_ds2);
                m.set(0, 0, 1, -dRga_dy*df_ds - Rga*d2f_dsdy);
                m.set(0, 1, 0, m.get(0, 0, 1));
                m.set(0, 1, 1, d2Rga_dy2 - d2Rga_dy2*f - 2.0*dRga_dy*df_dy - Rga*d2f_dy2);

                m.set(1, 0, 0, (Rob - Rgb)*d2f_ds2);
                m.set(1, 0, 1, (dRob_dy - dRgb_dy)*df_ds + (Rob - Rgb)*d2f_dsdy);
                m.set(1, 1, 0, m.get(1, 0, 1));
                m.set(1, 1, 1, d2Rgb_dy2 + (d2Rob_dy2 - d2Rgb_dy2)*f + 2.0*(dRob_dy - dRgb_dy)*df_dy + (Rob - Rgb)*d2f_dy2);

                m.set(2, 0, 0, Rod*d2f_ds2);
                m.set(2, 0, 1, dRod_dy*df_ds + Rod*d2f_dsdy);
                m.set(2, 1, 0, m.get(2, 0, 1));
                m.set(2, 1, 1, d2Rod_dy2*f + 2.0*dRod_dy*df_dy + Rod*d2f_dy2);
            }
        }
    }
    
    return degree;
}

int DeadVolatileVolatileGasFluxFunction::jet(const WaveState &state, JetMatrix &m, int degree) const {
    m.resize(3);

    double s = state(0);
    double y = state(1);
    double u = state(2);

    JetMatrix reduced;
    reduced_jet(state, reduced, degree);

    if (degree >= 0){
        // Copy the reduced jet and multiply it by u.
        //
        for (int i = 0; i < 3; i++) m.set(i, u*reduced.get(i));

        if (degree >= 1){
            // Copy the reduced jet's Jacobian and multiply it by u...
            //
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 2; j++){
                    m.set(i, j, u*reduced.get(i, j));
                }
            }

            // ...or  simply copy the reduced jet's function.
            //
            for (int i = 0; i < 3; i++) m.set(i, 2, reduced.get(i));

            if (degree >= 2){
                m.set(0, 0, 0, u*reduced.get(0, 0, 0));
                m.set(0, 0, 1, u*reduced.get(0, 0, 1));
                m.set(0, 0, 2, reduced.get(0, 0));
                
                m.set(0, 1, 0, m.get(0, 0, 1));
                m.set(0, 1, 1, u*reduced.get(0, 1, 1));
                m.set(0, 1, 2, reduced.get(0, 1));
                
                m.set(0, 2, 0, m.get(0, 0, 2));
                m.set(0, 2, 1, m.get(0, 1, 2));
                m.set(0, 2, 2, 0.0);
                
                // ====== //
                
                m.set(1, 0, 0, u*reduced.get(1, 0, 0));
                m.set(1, 0, 1, u*reduced.get(1, 0, 1));
                m.set(1, 0, 2, reduced.get(1, 0));
                
                m.set(1, 1, 0, m.get(1, 0, 1));
                m.set(1, 1, 1, u*reduced.get(1, 1, 1));
                m.set(1, 1, 2, reduced.get(1, 1));
                
                m.set(1, 2, 0, m.get(1, 0, 2));
                m.set(1, 2, 1, m.get(1, 1, 2));
                m.set(1, 2, 2, 0.0);
                
                // ====== //
                
                m.set(2, 0, 0, u*reduced.get(2, 0, 0));
                m.set(2, 0, 1, u*reduced.get(2, 0, 1));
                m.set(2, 0, 2, reduced.get(2, 0));
                
                m.set(2, 1, 0, m.get(2, 0, 1));
                m.set(2, 1, 1, u*reduced.get(2, 1, 1));
                m.set(2, 1, 2, reduced.get(2, 1));
                
                m.set(2, 2, 0, m.get(2, 0, 2));
                m.set(2, 2, 1, m.get(2, 1, 2));
                m.set(2, 2, 2, 0.0);
            }
        }
    }

    return degree;
}


#include "DeadVolatileVolatileGasAccumulationFunction.h"

DeadVolatileVolatileGasAccumulationFunction::DeadVolatileVolatileGasAccumulationFunction(Parameter *phi, DeadVolatileVolatileGasThermodynamics *th) : AccumulationFunction() {
    phi_parameter = phi;

    thermo = th;
}

DeadVolatileVolatileGasAccumulationFunction::~DeadVolatileVolatileGasAccumulationFunction(){
}


int DeadVolatileVolatileGasAccumulationFunction::reduced_jet(const WaveState &state, JetMatrix &m, int degree) const {
    m.resize(2, 3);

    double s = state(0);
    double y = state(1);

    // Get the updated parameters.
    //
    double phi = phi_parameter->value();

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
    
    if (degree >= 0){
        double Rga = Rga_jet.get(0);
        double Rgb = Rgb_jet.get(0);
        double Rob = Rob_jet.get(0);
        double Rod = Rod_jet.get(0);
        
        m.set(0, phi*Rga*(1.0 - s));
        m.set(1, phi*(Rgb + (Rob - Rgb)*s));
        m.set(2, phi*Rod*s);

        if (degree >= 1){
            double dRga_dy = Rga_jet.get(0, 0);
            double dRgb_dy = Rgb_jet.get(0, 0);
            
            double dRob_dy = Rob_jet.get(0, 0);
            double dRod_dy = Rod_jet.get(0, 0);
            
            m.set(0, 0, -phi*Rga);
            m.set(0, 1, phi*dRga_dy*(1.0 - s));

            m.set(1, 0, phi*(Rob - Rgb));
            m.set(1, 1, phi*(dRob_dy*s + dRgb_dy*(1.0 - s)));

            m.set(2, 0, phi*Rod);
            m.set(2, 1, phi*dRod_dy*s);

            if (degree >= 2){
                double d2Rga_dy2 = Rga_jet.get(0, 0, 0);
                double d2Rgb_dy2 = Rgb_jet.get(0, 0, 0);
                double d2Rob_dy2 = Rob_jet.get(0, 0, 0);
                double d2Rod_dy2 = Rod_jet.get(0, 0, 0);

                m.set(0, 0, 0, 0.0);
                m.set(0, 0, 1, -phi*dRga_dy);
                m.set(0, 1, 0, m.get(0, 0, 1));
                m.set(0, 1, 1, phi*d2Rga_dy2*(1.0 - s));

                m.set(1, 0, 0, 0.0);
                m.set(1, 0, 1, phi*(dRob_dy - dRgb_dy));
                m.set(1, 1, 0, m.get(1, 0, 1));
                m.set(1, 1, 1, phi*(d2Rgb_dy2 + (d2Rob_dy2 - d2Rgb_dy2)*s));

                m.set(2, 0, 0, 0.0);
                m.set(2, 0, 1, phi*dRod_dy);
                m.set(2, 1, 0, m.get(2, 0, 1));
                m.set(2, 1, 1, phi*d2Rod_dy2*s);
            }
        }
    }

    return degree;
}

int DeadVolatileVolatileGasAccumulationFunction::jet(const WaveState &state, JetMatrix &m, int degree) const {
    m.resize(3);
    
    JetMatrix reduced;
    reduced_jet(state, reduced, degree);
    
    if (degree >= 0){
        // Copy the reduced jet.
        //
        for (int i = 0; i < 3; i++) m.set(i, reduced.get(i));
        
        if (degree >= 1){
            m.set(0, 0, reduced.get(0, 0));
            m.set(0, 1, reduced.get(0, 1));
            m.set(0, 2, 0.0);

            m.set(1, 0, reduced.get(1, 0));
            m.set(1, 1, reduced.get(1, 1));
            m.set(1, 2, 0.0);

            m.set(2, 0, reduced.get(2, 0));
            m.set(2, 1, reduced.get(2, 1));
            m.set(2, 2, 0.0);
            
            if (degree >= 2){
                m.set(0, 0, 0, reduced.get(0, 0, 0));
                m.set(0, 0, 1, reduced.get(0, 0, 1));
                m.set(0, 0, 2, 0.0);

                m.set(0, 1, 0, reduced.get(0, 1, 0));
                m.set(0, 1, 1, reduced.get(0, 1, 1));
                m.set(0, 1, 2, 0.0);

                m.set(0, 2, 0, 0.0);
                m.set(0, 2, 1, 0.0);
                m.set(0, 2, 2, 0.0);

                //=====================//

                m.set(1, 0, 0, reduced.get(1, 0, 0));
                m.set(1, 0, 1, reduced.get(1, 0, 1));
                m.set(1, 0, 2, 0.0);

                m.set(1, 1, 0, reduced.get(1, 1, 0));
                m.set(1, 1, 1, reduced.get(1, 1, 1));
                m.set(1, 1, 2, 0.0);

                m.set(1, 2, 0, 0.0);
                m.set(1, 2, 1, 0.0);
                m.set(1, 2, 2, 0.0);

                //=====================//

                m.set(2, 0, 0, reduced.get(2, 0, 0));
                m.set(2, 0, 1, reduced.get(2, 0, 1));
                m.set(2, 0, 2, 0.0);

                m.set(2, 1, 0, reduced.get(2, 1, 0));
                m.set(2, 1, 1, reduced.get(2, 1, 1));
                m.set(2, 1, 2, 0.0);

                m.set(2, 2, 0, 0.0);
                m.set(2, 2, 1, 0.0);
                m.set(2, 2, 2, 0.0);
            }
            
            
        }
    }
    
    return degree;
}

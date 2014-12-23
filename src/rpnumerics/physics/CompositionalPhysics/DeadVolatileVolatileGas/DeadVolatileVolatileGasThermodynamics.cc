#include "DeadVolatileVolatileGasThermodynamics.h"

DeadVolatileVolatileGasThermodynamics::DeadVolatileVolatileGasThermodynamics(Parameter *B, Parameter *D, Parameter *mu_oB, Parameter *mu_oD, Parameter *mu_G, Parameter *rg, Parameter *re){
    B_parameter = B;
    D_parameter = D;

    mu_oB_parameter = mu_oB;
    mu_oD_parameter = mu_oD;
    mu_G_parameter  = mu_G;
    
    rg_parameter = rg;
    re_parameter = re;

    parameters_.push_back(B_parameter);
    parameters_.push_back(D_parameter);
    parameters_.push_back(mu_oB_parameter);
    parameters_.push_back(mu_oD_parameter);
    parameters_.push_back(mu_G_parameter);
    parameters_.push_back(rg_parameter);
    parameters_.push_back(re_parameter);

    info_auxiliaryfunction_ = std::string("Thermodynamics");
}

DeadVolatileVolatileGasThermodynamics::~DeadVolatileVolatileGasThermodynamics(){
}

void DeadVolatileVolatileGasThermodynamics::viscosity_ratio(int degree, double y, JetMatrix &R_jet){
    R_jet.resize(1);

    if (degree >= 0){
        double mu_oB = mu_oB_parameter->value();
        double mu_oD = mu_oD_parameter->value();
        double mu_G  = mu_G_parameter->value();

        double mB = mu_oB/mu_G;
        double mD = mu_oD/mu_G;

        double R = mB + (mD - mB)*y;
        R_jet.set(0, R);

        if (degree >= 1){
            double dR_dy = mD - mB;
            R_jet.set(0, 0, dR_dy);

            if (degree >= 2){
                double d2R_dy2 = 0.0;
                R_jet.set(0, 0, 0, d2R_dy2);
            }
        }
    }

    return;
}

void DeadVolatileVolatileGasThermodynamics::gas_molar_density_a(int degree, double y, JetMatrix &Rga_jet){
    Rga_jet.resize(1);
    
    if (degree >= 0){
        double rg = rg_parameter->value();
        double re = re_parameter->value();

        double Rga = rg - re*(1.0 - y);
        Rga_jet.set(0, Rga);

        if (degree >= 1){
            double dRga_dy = re;
            Rga_jet.set(0, 0, dRga_dy);

            if (degree >= 2){
                double d2Rga_dy2 = 0.0;
                Rga_jet.set(0, 0, 0, d2Rga_dy2);
            }
        }
    }
    
    return;
}

void DeadVolatileVolatileGasThermodynamics::gas_molar_density_b(int degree, double y, JetMatrix &Rgb_jet){
    Rgb_jet.resize(1);
    
    if (degree >= 0){
        double re = re_parameter->value();

        double Rgb = re*(1.0 - y);
        Rgb_jet.set(0, Rgb);

        if (degree >= 1){
            double dRgb_dy = -re;
            Rgb_jet.set(0, 0, dRgb_dy);

            if (degree >= 2){
                double d2Rgb_dy2 = 0.0;
                Rgb_jet.set(0, 0, 0, d2Rgb_dy2);
            }
        }
    }
    
    return;
}

void DeadVolatileVolatileGasThermodynamics::oil_molar_density(int degree, double y, JetMatrix &Ro_jet){
    Ro_jet.resize(1);

    if (degree >= 0){
        double B = B_parameter->value();
        double D = D_parameter->value();

        double L = D + (B - D)*y;
        double P = (D - B)/L;

        double Ro = B*D/L; // Oil density proper.
        Ro_jet.set(0, Ro);

        if (degree >= 1){
            double dRo_dy = P*Ro;
            Ro_jet.set(0, 0, dRo_dy);

            if (degree >= 2){
                double d2Ro_dy2 = 2.0*P*dRo_dy;
                Ro_jet.set(0, 0, 0, d2Ro_dy2);
            }
        }
    }

    return;
}

void DeadVolatileVolatileGasThermodynamics::oil_molar_density_b(int degree, const JetMatrix &Ro_jet, double y,  JetMatrix &Rob_jet){
    Rob_jet.resize(1);
    
    if (degree >= 0){
        double Ro = Ro_jet.get(0);
        double Rob = Ro*(1.0 - y);
        Rob_jet.set(0, Rob);

        if (degree >= 1){
            double dRo_dy = Ro_jet.get(0, 0);
            double dRob_dy = -Ro + (1.0 - y)*dRo_dy;
            Rob_jet.set(0, 0, dRob_dy);

            if (degree >= 2){
                double d2Ro_dy2 = Ro_jet.get(0, 0, 0);
                double d2Rob_dy2 = -2.0*dRo_dy + (1.0 - y)*d2Ro_dy2;
                Rob_jet.set(0, 0, 0, d2Rob_dy2);
            }
        }
    }
    
    return;
}

void DeadVolatileVolatileGasThermodynamics::oil_molar_density_d(int degree, const JetMatrix &Ro_jet, double y,  JetMatrix &Rod_jet){
    Rod_jet.resize(1);
    
    if (degree >= 0){
        double Ro = Ro_jet.get(0);
        double Rod = y*Ro;
        Rod_jet.set(0, Rod);

        if (degree >= 1){
            double dRo_dy = Ro_jet.get(0, 0);
            double dRod_dy = Ro + y*dRo_dy;
            Rod_jet.set(0, 0, dRod_dy);

            if (degree >= 2){
                double d2Ro_dy2 = Ro_jet.get(0, 0, 0);
                double d2Rod_dy2 = 2.0*dRo_dy + y*d2Ro_dy2;
                Rod_jet.set(0, 0, 0, d2Rod_dy2);
            }
        }
    }
    
    return;
}


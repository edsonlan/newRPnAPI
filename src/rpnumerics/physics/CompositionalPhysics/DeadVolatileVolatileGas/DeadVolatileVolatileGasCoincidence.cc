#include "DeadVolatileVolatileGasCoincidence.h"

DeadVolatileVolatileGasCoincidence::DeadVolatileVolatileGasCoincidence(DeadVolatileVolatileGasThermodynamics *th, 
                                                       DeadVolatileVolatileGasHydrodynamics *hy, 
                                                       Parameter *re, Parameter *rg, 
                                                       Parameter *phi) : Coincidence() {
    thermo = th;
    hydro  = hy;

    re_parameter = re;
    rg_parameter = rg;

    phi_parameter = phi;
}

DeadVolatileVolatileGasCoincidence::~DeadVolatileVolatileGasCoincidence(){
}

void DeadVolatileVolatileGasCoincidence::lambdas(const RealVector &p, double &lambda_s, double &lambda_e, double &lambda_diff) const {
    double s = p(0);
    double y = p(1);
    double u = p(2);

    double re  = re_parameter->value();
    double rg  = rg_parameter->value();
    double phi = phi_parameter->value();

    JetMatrix f_jet;
    hydro->fractional_flow(1, s, y, f_jet);
    double f     = f_jet.get(0);
    double df_ds = f_jet.get(0, 0);

    JetMatrix oildensity;
    thermo->oil_molar_density(0, y, oildensity); // Was: oil_density
    double Ro = oildensity.get(0);

//    double Gamma = rg*re*y*Ro;
//    double Psi   = Gamma - (rg - re*(1.0 - y))*(Ro*Ro);

    double Gamma = rg*re*y;
    double Psi   = Gamma - (rg - re*(1.0 - y))*Ro;

    double Omega = Gamma/Psi; 

    lambda_e = (u/phi)*(f - Omega)/(s - Omega);// std::cout << "Omega = " << Omega << std::endl;
    lambda_s = (u/phi)*df_ds;

    lambda_diff = lambda_s - lambda_e;

    return;
}

double DeadVolatileVolatileGasCoincidence::lambda_s(const RealVector &p) const {
    double ls, le, d;
    lambdas(p, ls, le, d);

    return ls;
}

double DeadVolatileVolatileGasCoincidence::lambda_e(const RealVector &p) const {
    double ls, le, d;
    lambdas(p, ls, le, d);

    return le;
}

double DeadVolatileVolatileGasCoincidence::lambda_diff(const RealVector &p) const {
    double ls, le, d;
    lambdas(p, ls, le, d);

    return d;
}


bool DeadVolatileVolatileGasCoincidence::extension_basis(const RealVector &u, double &fe, double &se) const {
    return true;
}


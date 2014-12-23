#include "DeadVolatileVolatileGasEvaporationExtension.h"

DeadVolatileVolatileGasEvaporationExtension::DeadVolatileVolatileGasEvaporationExtension(const DeadVolatileVolatileGasFluxFunction *f, 
                                                                         const DeadVolatileVolatileGasAccumulationFunction *a, 
                                                                         DeadVolatileVolatileGasCoincidence *g, Parameter *phi) : Extension(){
    flux = f;
    accumulation = a;
    god = g;
    phi_parameter = phi;
}

DeadVolatileVolatileGasEvaporationExtension::~DeadVolatileVolatileGasEvaporationExtension(){
}

int DeadVolatileVolatileGasEvaporationExtension::extension(const RealVector &p,  RealVector &ext_p){
    double sm = p(0);
    double ym = p(1);
    double um = p(2);

    double phi = phi_parameter->value();
    double lem = (god->lambda_e(p))*(phi/um);

    DeadVolatileVolatileGasThermodynamics* thermo = god->thermodynamics();
    JetMatrix R_jet;
    thermo->viscosity_ratio(0, ym, R_jet);
    double rm = R_jet.get(0);

    double A_Bhaskara = lem*(1.0 + rm)*((1.0 + rm)*sm*sm - 2.0*rm*sm + rm);
    double C_Bhaskara = rm*(lem*(1.0 + rm)*sm*sm - (1.0 + 2.0*rm*lem)*sm + rm*lem);
    double B_Bhaskara = -2.0*C_Bhaskara - rm;

    double x1 = 0.0, x2 = 0.0;
    int info = Utilities::Bhaskara(B_Bhaskara/A_Bhaskara, C_Bhaskara/A_Bhaskara, x1, x2);

//    std::cout << "A = " << A_Bhaskara << ", B = " << B_Bhaskara << ", C = " << C_Bhaskara << ", x1 = " << x1 << ", x2 = " << x2 << std::endl;
//    std::cout << "    info = " << info << std::endl;

    if (info == BHASKARA_COMPLEX_ROOTS) return EXTENSION_ERROR;

    if ((x1 < 0.0 || x1 > 1.0) && (x2 < 0.0 || x2 > 1.0)) return EXTENSION_ERROR;

    double root;
    if      ((x1 < 0.0 || x1 > 1.0) && (x2 >= 0.0 || x2 <= 1.0)) root = x2;
    else if ((x2 < 0.0 || x2 > 1.0) && (x1 >= 0.0 || x1 <= 1.0)) root = x1;
//    else root = (std::fabs(s_sigma_0 - x1) < std::fabs(s_sigma_0 - x2)) ? x1 : x2;
    else return EXTENSION_ERROR;

    ext_p.resize(3);

    // Extension of s_: s+.
    //
    ext_p(0) = root;

    // Extension of y_: y+ = y-.
    //
    ext_p(1) = ym;

    // Extension of u_: u+.
    // Two different points are needed now: + and -.
    //
    RealVector pp(0, 2, ext_p);
    RealVector pm(0, 2, p);

    JetMatrix Fpj, Fmj;
    flux->reduced_jet(pp, Fpj, 0);
    flux->reduced_jet(pm, Fmj, 0);
    RealVector Fp(Fpj.function());
    RealVector Fm(Fmj.function());

    JetMatrix Gpj, Gmj;
    accumulation->reduced_jet(pp, Gpj, 0);
    accumulation->reduced_jet(pm, Gmj, 0);
    RealVector Gp(Gpj.function());
    RealVector Gm(Gmj.function());

    int pos = 0;
    double maxF = 0.0;
    for (int i = 0; i < p.size(); i++){
        if (std::abs(Fp(i)) > maxF){
           pos = i;
           maxF = Fp(i);
        }
    }

    ext_p(2) = (lem*(Gp(pos) - Gm(pos)) + um*Fm(pos))/Fp(pos);

    return EXTENSION_OK;
}


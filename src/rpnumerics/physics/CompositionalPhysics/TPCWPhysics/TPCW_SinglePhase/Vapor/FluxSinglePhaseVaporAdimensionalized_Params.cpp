#include "FluxSinglePhaseVaporAdimensionalized_Params.h"

FluxSinglePhaseVaporAdimensionalized_Params::FluxSinglePhaseVaporAdimensionalized_Params(double mc, double mw, double Press, 
                                                                                         const char *hsigmaC_name,
                                                                                         double Tr, double Tw,
                                                                                         double Cr,
                                                                                         double cw,
                                                                                         double T_typical,
                                                                                         double Rho_typical,
                                                                                         double U_typical,
                                                                                         double h_typical) : FluxParams(11){
    component(0,  mc);
    component(1,  mw);
    component(2,  Press);
    component(3,  Tr);
    component(4,  Tw);
    component(5,  Cr);
    component(6,  cw);
    component(7,  T_typical);
    component(8,  Rho_typical);
    component(9,  U_typical);
    component(10, h_typical);

    hsigmaC = std::string(hsigmaC_name);
}

FluxSinglePhaseVaporAdimensionalized_Params::~FluxSinglePhaseVaporAdimensionalized_Params(){
}

const char *FluxSinglePhaseVaporAdimensionalized_Params::get_name(void){
    return hsigmaC.c_str();
}


#include "FoamPermeability.h"

FoamPermeability::FoamPermeability(Parameter *cnw, Parameter *cno, Parameter *cng,
                                   Parameter *nw,  Parameter *no,  Parameter *ng,
                                   ThreePhaseFlowSubPhysics *s) :  ThreePhaseFlowPermeability(s) {
    cnw_parameter = cnw;
    cno_parameter = cno;
    cng_parameter = cng;

    nw_parameter = nw;
    no_parameter = no;
    ng_parameter = ng;

    parameters_.push_back(cnw_parameter);
    parameters_.push_back(cno_parameter);
    parameters_.push_back(cng_parameter);

    info_auxiliaryfunction_ = std::string("Permeability");
}

FoamPermeability::~FoamPermeability(){
}

int FoamPermeability::PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water){
    double reduced_water;

    return PermeabilityWater_jet(state, degree, water, reduced_water);
}

int FoamPermeability::PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water, double &reduced_water){
    water.resize(2, 1);

    // Kludge.
    //
    double eps = 1e-6;

    double sw = state.component(0)/* + eps*/; // kludge to avoid problems with pow.

    double cnw = cnw_parameter->value();
    double cno = cno_parameter->value();
    double cng = cng_parameter->value();

    double inv_cn = 1.0/(1.0 - cnw - cno - cng);

    double swcnw = sw - cnw;

    if (swcnw <= 0.){ 
        if (degree >= 0){
            double kw = 0.;
            water.set(0, kw);

            // Reduced permeability.
            //
            reduced_water = 0.0;

            if (degree >= 1){
                double dkw_dsw = 0.;
                double dkw_dso = 0.;

                water.set(0, 0, dkw_dsw);
                water.set(0, 1, dkw_dso);

                if (degree >= 2){
                    double d2kw_dsw2  = 0.;
                    double d2kw_dswso = 0.;
                    double d2kw_dsosw = d2kw_dswso;
                    double d2kw_dso2  = 0.;

                    water.set(0, 0, 0, d2kw_dsw2);
                    water.set(0, 0, 1, d2kw_dswso);
                    water.set(0, 1, 0, d2kw_dsosw);
                    water.set(0, 1, 1, d2kw_dso2);
                }
            }
        }
    }
    else {
        if (degree >= 0){
            double scaled_sw = swcnw*inv_cn;
            double nw = nw_parameter->value();

            double pow2_kw = pow(scaled_sw, nw - 2.0);
            double pow1_kw = pow2_kw*scaled_sw;
            double kw      = pow1_kw*scaled_sw;

//            double scaled_sw = sw;
//            double nw = nw_parameter->value();

//            double pow2_kw = 1.0;
//            double pow1_kw = sw;
//            double kw      = sw*sw;

            water.set(0, kw);

            // Reduced permeability.
            //
            reduced_water = pow1_kw*inv_cn; // reduced_water = kw/swcnw; *** Aqui paramos. 27/10/2014

            if (degree >= 1){
                double dkw_dsw = nw*pow1_kw*inv_cn; // Rodrigo
                double dkw_dso = 0.;

                water.set(0, 0, dkw_dsw);
                water.set(0, 1, dkw_dso);

                if (degree >= 2){
                    double d2kw_dsw2  = nw*(nw - 1.0)*pow2_kw*inv_cn*inv_cn; // Rodrigo
                    double d2kw_dswso = 0.0;
                    double d2kw_dsosw = 0.0;
                    double d2kw_dso2  = 0.0;

                    water.set(0, 0, 0, d2kw_dsw2);
                    water.set(0, 0, 1, d2kw_dswso);
                    water.set(0, 1, 0, d2kw_dsosw);
                    water.set(0, 1, 1, d2kw_dso2);
                }
            }
        }
    }

    return degree;
}

int FoamPermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas){
    double reduced_gas;

    return PermeabilityGas_jet(state, degree, gas, reduced_gas);
}

int FoamPermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas, double &reduced_gas){
    gas.resize(2, 1);

    // Kludge.
    //
    double eps = 1e-6;

//    double sg =  1.0 - (state.component(0) + eps + state.component(1) + eps); // kludge to avoid problems with pow.
    double sg =  1.0 - (state.component(0) + state.component(1)); // kludge to avoid problems with pow.

    double cnw = cnw_parameter->value();
    double cno = cno_parameter->value();
    double cng = cng_parameter->value();

    double inv_cn = 1.0/(1.0 - cnw - cno - cng);

    double sgcng = sg - cng;

    if (sgcng <= 0.){ 
        if (degree >= 0){
            double kg = 0.;
            gas.set(0, kg);

            // Reduced permeability.
            //
            reduced_gas = 0.0;

            if (degree >= 1){
                double dkg_dsw = 0.;
                double dkg_dso = 0.;

                gas.set(0, 0, dkg_dsw);
                gas.set(0, 1, dkg_dso);

                if (degree >= 2){
                    double d2kg_dsw2  = 0.;
                    double d2kg_dswso = 0.;
                    double d2kg_dsosw = d2kg_dswso;
                    double d2kg_dso2  = 0.;

                    gas.set(0, 0, 0, d2kg_dsw2);
                    gas.set(0, 0, 1, d2kg_dswso);
                    gas.set(0, 1, 0, d2kg_dsosw);
                    gas.set(0, 1, 1, d2kg_dso2);
                }
            }
        }
    }
    else {
        if (degree >= 0){
            double scaled_sg = sgcng*inv_cn;
            double ng = nw_parameter->value();

            double pow2_kg = pow(scaled_sg, ng - 2.0);
            double pow1_kg = pow2_kg*scaled_sg;
            double kg      = pow1_kg*scaled_sg;

//            double scaled_sg = sg;
//            double ng = nw_parameter->value();

//            double pow2_kg = 1.0;
//            double pow1_kg = sg;
//            double kg      = sg*sg;

            gas.set(0, kg);

            // Reduced permeability.
            //
            reduced_gas = pow1_kg*inv_cn; // reduced_water = kw/swcnw; *** Aqui paramos. 27/10/2014

            if (degree >= 1){
                double dkg_dsw = -ng*pow1_kg*inv_cn; // Rodrigo
                double dkg_dso = dkg_dsw;

                gas.set(0, 0, dkg_dsw);
                gas.set(0, 1, dkg_dso);

                if (degree >= 2){
                    double d2kg_dsw2  = ng*(ng - 1.0)*pow2_kg*inv_cn*inv_cn; // Rodrigo
                    double d2kg_dswso = d2kg_dsw2;
                    double d2kg_dsosw = d2kg_dsw2;
                    double d2kg_dso2  = d2kg_dsw2;

                    gas.set(0, 0, 0, d2kg_dsw2);
                    gas.set(0, 0, 1, d2kg_dswso);
                    gas.set(0, 1, 0, d2kg_dsosw);
                    gas.set(0, 1, 1, d2kg_dso2);
                }
            }
        }
    }

    return degree;
}

//int FoamPermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas, double &reduced_gas){
//    gas.resize(2, 1);

//    // Kludge.
//    //
//    double eps = 1e-6;

////    double sg =  1.0 - (state.component(0) + eps + state.component(1) + eps); // kludge to avoid problems with pow.
//    double sg =  1.0 - (state.component(0) + state.component(1)); // kludge to avoid problems with pow.

//    double cnw = cnw_parameter->value();
//    double cno = cno_parameter->value();
//    double cng = cng_parameter->value();

//    double inv_cn = 1.0/(1.0 - cnw - cno - cng);

//    double sgcng = sg - cng;

//    if (sgcng <= 0.){ 
//        if (degree >= 0){
//            double kg = 0.;
//            gas.set(0, kg);

//            // Reduced permeability.
//            //
//            reduced_gas = 0.0;

//            if (degree >= 1){
//                double dkg_dsw = 0.;
//                double dkg_dso = 0.;

//                gas.set(0, 0, dkg_dsw);
//                gas.set(0, 1, dkg_dso);

//                if (degree >= 2){
//                    double d2kg_dsw2  = 0.;
//                    double d2kg_dswso = 0.;
//                    double d2kg_dsosw = d2kg_dswso;
//                    double d2kg_dso2  = 0.;

//                    gas.set(0, 0, 0, d2kg_dsw2);
//                    gas.set(0, 0, 1, d2kg_dswso);
//                    gas.set(0, 1, 0, d2kg_dsosw);
//                    gas.set(0, 1, 1, d2kg_dso2);
//                }
//            }
//        }
//    }
//    else {
//        if (degree >= 0){
//            double scaled_sg = sgcng*inv_cn;
//            double ng = nw_parameter->value();

//            double pow2_kg = pow(scaled_sg, ng - 2.0);
//            double pow1_kg = pow2_kg*scaled_sg;
//            double kg      = pow1_kg*scaled_sg;

////            double scaled_sg = sg;
////            double ng = nw_parameter->value();

////            double pow2_kg = 1.0;
////            double pow1_kg = sg;
////            double kg      = sg*sg;

//            gas.set(0, 1.0);

//            // Reduced permeability.
//            //
//            reduced_gas = pow1_kg*inv_cn; // reduced_water = kw/swcnw; *** Aqui paramos. 27/10/2014

//            if (degree >= 1){
//                double dkg_dsw = -ng*pow1_kg*inv_cn; // Rodrigo
//                double dkg_dso = dkg_dsw;

//                gas.set(0, 0, 0.0);
//                gas.set(0, 1, 0.0);

//                if (degree >= 2){
//                    double d2kg_dsw2  = ng*(ng - 1.0)*pow2_kg*inv_cn*inv_cn; // Rodrigo
//                    double d2kg_dswso = d2kg_dsw2;
//                    double d2kg_dsosw = d2kg_dsw2;
//                    double d2kg_dso2  = d2kg_dsw2;

//                    gas.set(0, 0, 0, 0.0);
//                    gas.set(0, 0, 1, 0.0);
//                    gas.set(0, 1, 0, 0.0);
//                    gas.set(0, 1, 1, 0.0);
//                }
//            }
//        }
//    }

//    return degree;
//}

int FoamPermeability::PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil){
    double reduced_oil;

    return PermeabilityOil_jet(state, degree, oil, reduced_oil);
}

int FoamPermeability::PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil, double &reduced_oil){
    oil.resize(2, 1);

    // Kludge.
    //
    double eps = 1e-6;

    double so = state.component(1)/* + eps*/; // kludge to avoid problems with pow.

    double cnw = cnw_parameter->value();
    double cno = cno_parameter->value();
    double cng = cng_parameter->value();

    double inv_cn = 1.0/(1.0 - cnw - cno - cng);

    double socno = so - cno;

    if (socno <= 0.){ 
        if (degree >= 0){
            double ko = 0.;
            oil.set(0, ko);

            // Reduced permeability.
            //
            reduced_oil = 0.0;

            if (degree >= 1){
                double dko_dsw = 0.;
                double dko_dso = 0.;

                oil.set(0, 0, dko_dsw);
                oil.set(0, 1, dko_dso);

                if (degree >= 2){
                    double d2ko_dsw2  = 0.;
                    double d2ko_dswso = 0.;
                    double d2ko_dsosw = d2ko_dswso;
                    double d2ko_dso2  = 0.;

                    oil.set(0, 0, 0, d2ko_dsw2);
                    oil.set(0, 0, 1, d2ko_dswso);
                    oil.set(0, 1, 0, d2ko_dsosw);
                    oil.set(0, 1, 1, d2ko_dso2);
                }
            }
        }
    }
    else {
        if (degree >= 0){
            double scaled_so = socno*inv_cn;
            double no = no_parameter->value();

            double pow2_ko = pow(scaled_so, no - 2.0);
            double pow1_ko = pow2_ko*scaled_so;
            double ko      = pow1_ko*scaled_so;

//            double scaled_so = so;
//            double no = no_parameter->value();

//            double pow2_ko = 1.0;
//            double pow1_ko = so;
//            double ko      = so*so;

            oil.set(0, ko);

            // Reduced permeability.
            //
            reduced_oil = pow1_ko*inv_cn;

            if (degree >= 1){
                double dko_dsw = 0.0;
                double dko_dso = no*pow1_ko*inv_cn;

                oil.set(0, 0, dko_dsw);
                oil.set(0, 1, dko_dso);

                if (degree >= 2){
                    double d2ko_dsw2  = 0.0;
                    double d2ko_dswso = 0.0;
                    double d2ko_dsosw = 0.0;
                    double d2ko_dso2  = no*(no - 1.0)*pow2_ko*inv_cn*inv_cn;

                    oil.set(0, 0, 0, d2ko_dsw2);
                    oil.set(0, 0, 1, d2ko_dswso);
                    oil.set(0, 1, 0, d2ko_dsosw);
                    oil.set(0, 1, 1, d2ko_dso2);
                }
            }
        }
    }

    return degree;
}

void FoamPermeability::reduced_permeability(const RealVector &state, RealVector &reduced){
    double rw, ro, rg;

    JetMatrix water, oil, gas;
    PermeabilityWater_jet(state, 0, water, rw);
    PermeabilityGas_jet(state, 0, gas, rg);
    PermeabilityOil_jet(state, 0, oil, ro);

    reduced.resize(3);
    reduced(0) = rw;
    reduced(1) = ro;
    reduced(2) = rg;

    return;
}


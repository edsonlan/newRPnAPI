#include "Brooks_CoreyPermeability.h"

double Brooks_CoreyPermeability::power(double x, double y){
    double p = 1.0;

    for (int i = 0; i < (int)y; i++) p *= x;

    return p;
}

Brooks_CoreyPermeability::Brooks_CoreyPermeability(Parameter *lambda, Parameter *cnw, Parameter *cno, Parameter *cng, ThreePhaseFlowSubPhysics *s) :  ThreePhaseFlowPermeability(s) {
    lambda_parameter = lambda;
    cnw_parameter = cnw;
    cno_parameter = cno;
    cng_parameter = cng;

    parameters_.push_back(lambda_parameter);
    parameters_.push_back(cnw_parameter);
    parameters_.push_back(cno_parameter);
    parameters_.push_back(cng_parameter);

    info_auxiliaryfunction_ = std::string("Permeability");
}

Brooks_CoreyPermeability::~Brooks_CoreyPermeability(){
}

int Brooks_CoreyPermeability::PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water){
    double reduced_water;

    return PermeabilityWater_jet(state, degree, water, reduced_water);
}

int Brooks_CoreyPermeability::PermeabilityWater_jet(const RealVector &state, int degree, JetMatrix &water, double &reduced_water){
    // std::cout << "Water begins." << std::endl;

    water.resize(2, 1);

    // Kludge.
    //
    double eps = 1e-6;

    double sw = state.component(0) + eps; // kludge to avoid problems with pow.
    double so = state.component(1) + eps; // kludge to avoid problems with pow.
//    double sg = 1.0 - sw - so;

    double cnw = cnw_parameter->value();

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
                    double d2kw_dsw2 = 0.;
                    double d2kw_dswso = 0.;
                    double d2kw_dsosw = d2kw_dswso;
                    double d2kw_dso2 = 0.;

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
            double lambda_s = lambda_parameter->value();
            double a = (2.0 + lambda_s)/lambda_s;
            // // std::cout << "Water. a = " << a << ", swcnw = " << swcnw << std::endl;

//            double swcnwE2 = power(swcnw, a);
            double swcnwE2 = pow(swcnw, a);

            // // std::cout << "    After pow." << std::endl;

            double swcnwE1 = swcnw*swcnwE2;

            double kw = swcnw*swcnwE1;

            water.set(0, kw);

            // Reduced permeability.
            //
            reduced_water = swcnwE1;

            if (degree >= 1){
                double a2 = a + 2.0;
                double dkw_dsw = a2*swcnwE1;
                double dkw_dso = 0.;

                water.set(0, 0, dkw_dsw);
                water.set(0, 1, dkw_dso);

                if (degree >= 2){
                    double a1 = a + 1.0;

                    double d2kw_dsw2 = a2*a1*swcnwE2;
                    double d2kw_dswso = 0.;
                    double d2kw_dsosw = d2kw_dswso;
                    double d2kw_dso2 = 0.;

                    water.set(0, 0, 0, d2kw_dsw2);
                    water.set(0, 0, 1, d2kw_dswso);
                    water.set(0, 1, 0, d2kw_dsosw);
                    water.set(0, 1, 1, d2kw_dso2);
                }
            }
        }
    }

    // std::cout << "    Water ends." << std::endl;

    return degree;
}

int Brooks_CoreyPermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas){
    double reduced_gas;

    return PermeabilityGas_jet(state, degree, gas, reduced_gas);
}

int Brooks_CoreyPermeability::PermeabilityGas_jet(const RealVector &state, int degree, JetMatrix &gas, double &reduced_gas){
    // std::cout << "Gas begins." << std::endl;

    gas.resize(2, 1);


    // Kludge.
    //
    double eps = 1e-6;

    double sw = state.component(0) + eps; // kludge to avoid problems with pow.
    double so = state.component(1) + eps; // kludge to avoid problems with pow.

//    double sw = state.component(0);
//    double so = state.component(1);
    double sg = 1.0 - sw - so;

    double cng = cng_parameter->value();
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
                    double d2kg_dsw2 = 0.;
                    double d2kg_dswso = 0.;
                    double d2kg_dsosw = d2kg_dswso;
                    double d2kg_dso2 = 0.;

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
            double lambda_s = lambda_parameter->value();
            double a = (2.0 + lambda_s)/lambda_s;

            // // std::cout << "Gas. Before pow1: 1.0 - sgcng = " << 1.0 - sgcng << ", a - 2.0 = "<< a - 2.0 << std::endl;
//            double swsoE2 = power(1.0 - sgcng, a - 2.0);
            double swsoE2 = pow(1.0 - sgcng, a - 2.0);
            // // std::cout << "    After pow1." << std::endl;

            double swsoE1 = swsoE2*(1.0 - sgcng);
            double Term1 = 1.0 - (1.0 - sgcng)*swsoE1;

            double kg = sgcng*sgcng*Term1;
            gas.set(0, kg);

            // Reduced permeability.
            //
            reduced_gas = sgcng*Term1;

            if (degree >= 1){
                double dkg_dsw = -2.0*sgcng*Term1 - a*sgcng*sgcng*swsoE1;
                double dkg_dso = dkg_dsw;

                gas.set(0, 0, dkg_dsw);
                gas.set(0, 1, dkg_dso);

                if (degree >= 2){
                    double d2kg_dsw2  = 2.0*Term1 + 4.0*sgcng*a*swsoE1 - a*(a - 1.0)*sgcng*sgcng*swsoE2;
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

    // std::cout << "    Gas ends." << std::endl;

    return degree;
}

int Brooks_CoreyPermeability::PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil){
    double reduced_oil;

    return PermeabilityOil_jet(state, degree, oil, reduced_oil);
}

int Brooks_CoreyPermeability::PermeabilityOil_jet(const RealVector &state, int degree, JetMatrix &oil, double &reduced_oil){
    // std::cout << "Oil begins." << std::endl;

    oil.resize(2, 1);

    // Kludge.
    //
    double eps = 1e-6;

    double sw = state.component(0) + eps; // kludge to avoid problems with pow.
    double so = state.component(1) + eps; // kludge to avoid problems with pow.

//    double sw = state.component(0);
//    double so = state.component(1);
    double sg = 1.0 - sw - so;

    double cnw = cnw_parameter->value();
    double cno = cno_parameter->value();

    double swcnw = sw - cnw;
    double socno = so - cno;

    if (socno <= 0.){ 
        if (degree >= 0){
            double ko = 0.0;
            oil.set(0, ko);

            // Reduced permeability.
            //
            reduced_oil = 0.0;

            if (degree >= 1){
                double dko_dsw = 0.0;
                double dko_dso = 0.0;

                oil.set(0, 0, dko_dsw);
                oil.set(0, 1, dko_dso);

                if (degree >= 2){
                    double d2ko_dsw2 = 0.0;
                    double d2ko_dswso = 0.0;
                    double d2ko_dsosw = d2ko_dswso;
                    double d2ko_dso2 = 0.;

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
            double lambda_s = lambda_parameter->value();
            double a       = (2.0 + lambda_s)/lambda_s;

            // // std::cout << "Oil. Before pow, swcnw + socno = " << swcnw + socno << ", a - 1.0 = " << a - 1.0 << std::endl;
//            double swsoE1  = power(swcnw + socno, a - 1.0);
            double swsoE1  = pow(swcnw + socno, a - 1.0);
            // // std::cout << "    After pow." << std::endl;

            double swsoE0  = (swcnw + socno)*swsoE1;

            // // std::cout << "Oil, before 1.0 - power(swcnw, a), swcnw = " << swcnw << ", a = " << a << std::endl;
            double swa1    = 1.0 - power(swcnw, a);
            // // std::cout << "    After pow." << std::endl;

            double swsoE_1 = (swcnw + socno)*swsoE0;

            // // std::cout << "Oil, before power(swcnw, a - 2.0), swcnw = " << swcnw << ", a - 2.0 = " << a - 2.0 << std::endl;
            double pswE2 = power(swcnw, a - 2.0);
            // // std::cout << "    After pow." << std::endl;

            double pswE1 = swcnw*pswE2;

            double ko = socno*(1. - swcnw)*swa1*swsoE_1;
            oil.set(0, ko);

            // Reduced permeability.
            //
            reduced_oil = (1. - swcnw)*swa1*swsoE_1;

            if (degree >= 1){
                double a1 = a + 1.0;
                double a2 = a1 + 1.0;
		
                double dko_dsw = -socno*swa1*swsoE_1 - a*socno*pswE1*(1.0 - swcnw)*swsoE_1 + a1*socno*(1.0 - swcnw)*swa1*swsoE0;
                double dko_dso = (swcnw + a2*socno)*(1.0 - swcnw)*swa1*swsoE0;

                oil.set(0, 0, dko_dsw);
                oil.set(0, 1, dko_dso);

                if (degree >= 2){
                    double d2ko_dsw2  = socno*swsoE1*(
                                                      (swcnw + socno)*(
                                                                       2.0*a*pswE1*(swcnw + socno) - 
                                                                       2.0*a1*swa1 - 
                                                                       a*(a - 1.0)*pswE2*(1.0 - swcnw)*(swcnw + socno) - 
                                                                       2.0*a*a1*pswE1*(1.0 - swcnw)
                                                                      ) 
                                                      + a*a1*(1.0 - swcnw)*swa1
                                                     );
                    double d2ko_dswso = (1.0 - swcnw)*swa1*swsoE0 + (swcnw + a2*socno)*(
                                                                                        a*(1.0 - swcnw)*swa1*swsoE1 - 
                                                                                        swsoE0*(swa1 + a*pswE1*(1.0 - swcnw))
                                                                                       );
                    double d2ko_dsosw = d2ko_dswso;
                    double d2ko_dso2  = (1.0 - swcnw)*swa1*swsoE1*(a2*(swcnw + socno) + a*(swcnw + a2*socno));

                    oil.set(0, 0, 0, d2ko_dsw2);
                    oil.set(0, 0, 1, d2ko_dswso);
                    oil.set(0, 1, 0, d2ko_dsosw);
                    oil.set(0, 1, 1, d2ko_dso2);
                }
            }
        }
    }

    // std::cout << "    Oil ends." << std::endl;

    return degree;
}

void Brooks_CoreyPermeability::reduced_permeability(const RealVector &state, RealVector &reduced){
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


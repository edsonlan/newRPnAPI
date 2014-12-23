#include "ICDOWChemistry.h"

ICDOWChemistry::ICDOWChemistry(){
}

ICDOWChemistry::~ICDOWChemistry(){
}

// x = Hydrogen contents.
//
void ICDOWChemistry::aqueous_carbon_excess_oxygen(double x, int degree, JetMatrix &rho_CO_jet){
    rho_CO_jet.resize(1);

    if (degree >= 0){
        rho_CO_jet.set(0, 5.625319783841633E2/(x*3.30182732201303E1+sinh(x)+2.31145688182898E2));

        if (degree >= 1){
            rho_CO_jet.set(0, 0, (cosh(x)+3.30182732201303E1)*1.0/pow(x*3.30182732201303E1+sinh(x)+2.31145688182898E2,2.0)*(-5.625319783841633E2));

            if (degree >= 2){
                rho_CO_jet.set(0, 0, 0, sinh(x)*1.0/pow(x*3.30182732201303E1+sinh(x)+2.31145688182898E2,2.0)*(-5.625319783841633E2)+pow(cosh(x)+3.30182732201303E1,2.0)*1.0/pow(x*3.30182732201303E1+sinh(x)+2.31145688182898E2,3.0)*1.125063956768327E3
                              );
            }
        }
    }

    return;
}

// x = Hydrogen contents.
//
void ICDOWChemistry::aqueous_hydrogen(double x, int degree, JetMatrix &rho_H_jet){
    rho_H_jet.resize(1);

    if (degree >= 0){
        rho_H_jet.set(0, exp(1.0751E1/(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3)-1.0/1.0E2));

        if (degree >= 1){
            rho_H_jet.set(0, 0, exp(1.0751E1/(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3)-1.0/1.0E2)*(cosh(x*1.50509958419512)*1.50509958419512+sinh(x)*3.8642421361418E1)*1.0/pow(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3,2.0)*(-1.0751E1)
                         );

            if (degree >= 2){
                rho_H_jet.set(0, 0, 0, exp(1.0751E1/(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3)-1.0/1.0E2)*(sinh(x*1.50509958419512)*2.265324758344323+cosh(x)*3.8642421361418E1)*1.0/pow(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3,2.0)*(-1.0751E1)+exp(1.0751E1/(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3)-1.0/1.0E2)*pow(cosh(x*1.50509958419512)*1.50509958419512+sinh(x)*3.8642421361418E1,2.0)*1.0/pow(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3,3.0)*2.1502E1+exp(1.0751E1/(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3)-1.0/1.0E2)*pow(cosh(x*1.50509958419512)*1.50509958419512+sinh(x)*3.8642421361418E1,2.0)*1.0/pow(sinh(x*1.50509958419512)+cosh(x)*3.8642421361418E1-6.4299174481606E3,4.0)*1.15584001E2
                             );
            }
        }
    }

    return;
}

// x = Hydrogen contents.
//
void ICDOWChemistry::carbon_in_oil(double x, int degree, JetMatrix &rho_C_oil_jet){
    rho_C_oil_jet.resize(1);

    if (degree >= 0){
        rho_C_oil_jet.set(0, exp(tanh(sinh(x)*9.4728E-4)*(-3.29E2/5.0E1))*(1.0/7.0E2));

        if (degree >= 1){
            rho_C_oil_jet.set(0, 0, exp(tanh(sinh(x)*9.4728E-4)*(-3.29E2/5.0E1))*cosh(x)*(pow(tanh(sinh(x)*9.4728E-4),2.0)-1.0)*8.904432E-6
                             );

            if (degree >= 2){
                rho_C_oil_jet.set(0, 0, 0, exp(tanh(sinh(x)*9.4728E-4)*(-3.29E2/5.0E1))*pow(cosh(x),2.0)*pow(pow(tanh(sinh(x)*
                                           9.4728E-4),2.0)-1.0,2.0)*5.55022364698368E-8+exp(tanh(sinh(x)*9.4728E-4)*(-3.29E2/5.0E1))*sinh(x)*
                                           (pow(tanh(sinh(x)*9.4728E-4),2.0)-1.0)*8.904432E-6-exp(tanh(sinh(x)*9.4728E-4)*(-3.29E2/5.0E1))*
                                           tanh(sinh(x)*9.4728E-4)*pow(cosh(x),2.0)*(pow(tanh(sinh(x)*9.4728E-4),2.0)-1.0)*1.686998068992E-8
                                 );
            }
        }
    }

    return;
}


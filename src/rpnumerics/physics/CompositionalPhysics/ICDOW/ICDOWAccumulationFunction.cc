#include "ICDOWAccumulationFunction.h"

ICDOWAccumulationFunction::ICDOWAccumulationFunction(Parameter *phi, ICDOWChemistry *ch): AccumulationFunction()
                                                                                          
 {
    phi_parameter = phi;

    chemistry = ch;

    std::cout << "Accumulation = " << (void*)this << std::endl;
}

ICDOWAccumulationFunction::~ICDOWAccumulationFunction(){
}


int ICDOWAccumulationFunction::reduced_jet(const WaveState &state, JetMatrix &m, int degree) const {
    m.resize(2, 3);

    double sw  = state(0); // Water saturation.
    double hyd = state(1); // Hydrogen concentration.

    double so = 1.0 - sw;

    // Get the updated parameter.
    //
    double phi = phi_parameter->value();

    // rho_a,C-O.
    //
    JetMatrix rho_CO_jet;
    chemistry->aqueous_carbon_excess_oxygen(hyd, degree, rho_CO_jet);

    // rho_a,H(1).
    //
    JetMatrix rho_H_jet;
    chemistry->aqueous_hydrogen(hyd, degree, rho_H_jet);

    // rho carbon in oil.
    //
    JetMatrix rho_C_oil_jet;
    chemistry->carbon_in_oil(hyd, degree, rho_C_oil_jet);
    
    if (degree >= 0){
        m.set(0, phi*sw*rho_CO_jet.get(0));
        m.set(1, phi*sw*rho_H_jet.get(0));
        m.set(2, phi*so*rho_C_oil_jet.get(0));

        if (degree >= 1){
            m.set(0, 0, phi*rho_CO_jet.get(0));
            m.set(0, 1, phi*sw*rho_CO_jet.get(0, 0));

            m.set(1, 0, phi*rho_H_jet.get(0));
            m.set(1, 1, phi*sw*rho_H_jet.get(0, 0));

            m.set(2, 0, -phi*rho_C_oil_jet.get(0));
            m.set(2, 1, phi*so*rho_C_oil_jet.get(0, 0));

            if (degree >= 2){
                m.set(0, 0, 0, 0.0);
                m.set(0, 0, 1, phi*rho_CO_jet.get(0, 0));
                m.set(0, 1, 0, m.get(0, 0, 1));
                m.set(0, 1, 1, phi*sw*rho_CO_jet.get(0, 0, 0));

                m.set(1, 0, 0, 0.0);
                m.set(1, 0, 1, phi*rho_H_jet.get(0, 0));
                m.set(1, 1, 0, m.get(1, 0, 1));
                m.set(1, 1, 1, phi*sw*rho_H_jet.get(0, 0, 0));

                m.set(2, 0, 0, 0.0);
                m.set(2, 0, 1, -phi*rho_C_oil_jet.get(0, 0));
                m.set(2, 1, 0, m.get(2, 0, 1));
                m.set(2, 1, 1, phi*so*rho_C_oil_jet.get(0, 0, 0));
            }
        }
    }

    return degree;
}

int ICDOWAccumulationFunction::jet(const WaveState &state, JetMatrix &m, int degree) const {
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

#include "ICDOWFluxFunction.h"

ICDOWFluxFunction::ICDOWFluxFunction(ICDOWChemistry *ch, ICDOWHydrodynamics *hy): FluxFunction()
                                                                                  
{
    chemistry = ch;
    hydro     = hy;
}

ICDOWFluxFunction::~ICDOWFluxFunction(){
}

int ICDOWFluxFunction::reduced_jet(const WaveState &state, JetMatrix &m, int degree) const {
    m.resize(2, 3);

    double sw  = state(0);
    double hyd = state(1);

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

    // Hydrodynamics.
    //
    JetMatrix fw_jet;
    hydro->water_fractional_flow(sw, degree, fw_jet);

    JetMatrix fo_jet;
    hydro->oil_fractional_flow(sw, degree, fo_jet);
    
    if (degree >= 0){
        double fw = fw_jet.get(0);
        double fo = fo_jet.get(0);

//        m.set(0, fw/* *rho_CO_jet.get(0)*/);
        m.set(0, fw*rho_CO_jet.get(0));
        m.set(1, fw*rho_H_jet.get(0));
        m.set(2, fo*rho_C_oil_jet.get(0));
        
        if (degree >= 1){
            double dfw_dsw = fw_jet.get(0, 0);
            double dfo_dsw = fo_jet.get(0, 0);

//            m.set(0, 0, dfw_dsw/* *rho_CO_jet.get(0)*/);
//            m.set(0, 1, 0.0 /*fw*rho_CO_jet.get(0, 0)*/);
            m.set(0, 0, dfw_dsw*rho_CO_jet.get(0));
            m.set(0, 1, fw*rho_CO_jet.get(0, 0));


            m.set(1, 0, dfw_dsw*rho_H_jet.get(0));
            m.set(1, 1, fw*rho_H_jet.get(0, 0));

            m.set(2, 0, dfo_dsw*rho_C_oil_jet.get(0));
            m.set(2, 1, fo*rho_C_oil_jet.get(0, 0));
            
            if (degree >= 2){
                double d2fw_dsw2 = fw_jet.get(0, 0, 0);
                double d2fo_dsw2 = fo_jet.get(0, 0, 0);

                m.set(0, 0, 0, d2fw_dsw2*rho_CO_jet.get(0));
                m.set(0, 0, 1, dfw_dsw*rho_CO_jet.get(0, 0));
                m.set(0, 1, 0, m.get(0, 0, 1));
                m.set(0, 1, 1, fw*rho_CO_jet.get(0, 0, 0));

                m.set(1, 0, 0, d2fw_dsw2*rho_H_jet.get(0));
                m.set(1, 0, 1, dfw_dsw*rho_H_jet.get(0, 0));
                m.set(1, 1, 0, m.get(1, 0, 1));
                m.set(1, 1, 1, fw*rho_H_jet.get(0, 0, 0));

                m.set(2, 0, 0, d2fo_dsw2*rho_C_oil_jet.get(0));
                m.set(2, 0, 1, dfo_dsw*rho_C_oil_jet.get(0, 0));
                m.set(2, 1, 0, m.get(2, 0, 1));
                m.set(2, 1, 1, fo*rho_C_oil_jet.get(0, 0, 0));
            }
        }
    }
    
    return degree;
}

int ICDOWFluxFunction::jet(const WaveState &state, JetMatrix &m, int degree) const {
    m.resize(3);

    double sw  = state(0);
    double hyd = state(1);
    double u   = state(2);

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


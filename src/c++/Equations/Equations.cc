#include "Equations.h"

Equations::Equations(){
}

Equations::~Equations(){
}

int Equations::constrained_conservation_laws_jet(const RealVector &p, DoubleMatrix &DFcompleted, DoubleMatrix &DGcompleted){
    // Use the flux and the accumulation.
    //
    JetMatrix Fjet, Gjet, Cjet;

    int degree = 1;
    compute(p, degree, Fjet, Gjet, Cjet);

    DoubleMatrix Fjac = Fjet.Jacobian();
    DoubleMatrix Gjac = Gjet.Jacobian();

    // Here it is not said, but the matrices returned by this method are square.
    // In other words: Feqs + Ceqs == Fvars.
    //
    int rows = number_equations + number_constraints; // Merely mnemonic, since rows = number_variables.
    int cols = number_variables;

    DFcompleted.resize(rows, cols);
    DGcompleted.resize(rows, cols);

    for (int j = 0; j < number_variables; j++){
        for (int i = 0; i < number_equations; i++){
            DFcompleted(i, j) = Fjac(i, j);
            DGcompleted(i, j) = Gjac(i, j);
        }
    }

    if (number_constraints > 0){
        DoubleMatrix Cjac = Cjet.Jacobian();

        for (int j = 0; j < number_variables; j++){
            for (int i = 0; i < number_constraints; i++){
                DFcompleted(i + number_equations, j) = Cjac(i, j);
                DGcompleted(i + number_equations, j) = 0.0;
            }
        }
    }

    return CONSTRAINEDCONVECTION_OK;
}


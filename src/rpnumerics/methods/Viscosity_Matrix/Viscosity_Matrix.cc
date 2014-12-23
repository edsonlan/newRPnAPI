#include "Viscosity_Matrix.h"

void Viscosity_Matrix::fill_viscous_matrix(const RealVector &p, ViscosityJetMatrix &m) const {
    fill_viscous_matrix(p, m, 0);
    return;
}

void Viscosity_Matrix::fill_viscous_matrix(const RealVector &p, ViscosityJetMatrix &m, int degree) const {
    if (degree >= 0){
        m.M() = DoubleMatrix::eye(3); // Identity
//        m.M()(2, 2) = 0.0;

        if (degree >= 1){
            printf("Viscosity_Matrix::fill_viscous_matrix(): This method is not prepared yet to deal with degrees higher than 0!\n");
            exit(0);

            if (degree >= 2){
            }
        }
    }

    return;
}


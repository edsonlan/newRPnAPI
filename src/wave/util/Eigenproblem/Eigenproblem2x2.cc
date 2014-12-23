#include "Eigenproblem2x2.h"

Eigenproblem2x2::Eigenproblem2x2(): Eigenproblem(){
}

Eigenproblem2x2::~Eigenproblem2x2(){
}

// Standard problem.
//
void Eigenproblem2x2::find_eigenpair(const DoubleMatrix &A, int index, Eigenpair &ep){
}

void Eigenproblem2x2::find_eigenpairs(const DoubleMatrix &A, std::vector<Eigenpair> &eps){
}

void Eigenproblem2x2::find_eigenvalue(const DoubleMatrix &A, int index, Eigenvalue &ev){
    double a = 1.0;
    double b = -(A(0, 0) + A(1, 1));
    double c = A(0, 0)*A(1, 1) - A(1, 0)*A(0, 1);

    if (b >= 0.0){
    }
    else {
    }

    return;
}

void Eigenproblem2x2::find_eigenvalues(const DoubleMatrix &A, std::vector<Eigenvalue> &evs){
    evs.resize(2);

    double a = 1.0;
    double b = -(A(0, 0) + A(1, 1));
    double c = A(0, 0)*A(1, 1) - A(1, 0)*A(0, 1);

    double Delta = b*b - 4.0*a*c;
    double sqrtDelta = sqrt(std::abs(Delta));

    double inv2a = 1.0/(2.0*a);

    if (Delta > 0){
        for (int i = 0; i < 2; i++){
            evs[i].is_real   = true;
            evs[i].imaginary = 0.0;
        }

        if (b >= 0.0){
            evs[0].real = (-b - sqrtDelta)*inv2a;
            evs[1].real = -2.0*c/(b + sqrtDelta);
        }
        else {
            evs[0].real = 2.0*c/(-b + sqrtDelta);
            evs[1].real = (-b + sqrtDelta)*inv2a;
        }

    } 
    else {
        evs[0].is_real = evs[1].is_real = false;

        evs[0].real = evs[1].real = -b*inv2a;

        evs[1].imaginary = sqrtDelta*inv2a;
        evs[0].imaginary = -evs[1].imaginary;
    }

    return;
}

// Generalized problem.
//
void Eigenproblem2x2::find_eigenpair(const DoubleMatrix &A, const DoubleMatrix &B, int index, Eigenpair &ep){
}

void Eigenproblem2x2::find_eigenpairs(const DoubleMatrix &A, const DoubleMatrix &B, std::vector<Eigenpair> &eps){
}

void Eigenproblem2x2::find_eigenvalue(const DoubleMatrix &A, const DoubleMatrix &B, int index, Eigenvalue &ev){
}

void Eigenproblem2x2::find_eigenvalues(const DoubleMatrix &A, const DoubleMatrix &B, std::vector<Eigenvalue> &evs){
}


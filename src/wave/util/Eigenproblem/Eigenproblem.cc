#include "Eigenproblem.h"

Eigenproblem::Eigenproblem(){
}

Eigenproblem::~Eigenproblem(){
}

// Standard problem.
//
void Eigenproblem::find_eigenpair(const DoubleMatrix &A, int index, Eigenpair &ep){
}
void Eigenproblem::find_eigenpairs(const DoubleMatrix &A, std::vector<Eigenpair> &eps){
}

void Eigenproblem::find_eigenvalue(const DoubleMatrix &A, int index, Eigenvalue &ev){
}
void Eigenproblem::find_eigenvalues(const DoubleMatrix &A, std::vector<Eigenvalue> &evs){
}

// Generalized problem.
//
void Eigenproblem::find_eigenpair(const DoubleMatrix &A, const DoubleMatrix &B, int index, Eigenpair &ep){
}
void Eigenproblem::find_eigenpairs(const DoubleMatrix &A, const DoubleMatrix &B, std::vector<Eigenpair> &eps){
}

void Eigenproblem::find_eigenvalue(const DoubleMatrix &A, const DoubleMatrix &B, int index, Eigenvalue &ev){
}
void Eigenproblem::find_eigenvalues(const DoubleMatrix &A, const DoubleMatrix &B, std::vector<Eigenvalue> &evs){
}


#include "Eigenpair.h"

void Eigenpair::copy(const Eigenpair &orig){
    is_real = orig.is_real;

    eigenvalue = orig.eigenvalue;
    right_eigenvector = orig.right_eigenvector;
    left_eigenvector = orig.left_eigenvector;

    return;
}

Eigenpair::Eigenpair(){
}

Eigenpair::Eigenpair(int n){
    right_eigenvector.real.resize(n);
    right_eigenvector.imaginary.resize(n);

    left_eigenvector.real.resize(n);
    left_eigenvector.imaginary.resize(n);
}

Eigenpair::Eigenpair(const Eigenpair &orig){
    copy(orig);
}

Eigenpair::Eigenpair(const Eigenpair *orig){
    copy(*orig);
}

Eigenpair::~Eigenpair(){
}


Eigenpair Eigenpair::operator=(const Eigenpair &orig){
    if (&orig != this) copy(orig);

    return *this;
}


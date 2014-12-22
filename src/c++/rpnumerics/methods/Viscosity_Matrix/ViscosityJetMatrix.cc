#include "ViscosityJetMatrix.h"
#include <iostream>

void ViscosityJetMatrix::init(){
    M_.resize(nrows_, ncols_);
    //JM_.resize(nrows_, ncols_); // Todo: reserve space with for cycle

    return;
}

ViscosityJetMatrix::ViscosityJetMatrix() {
    nrows_ = ncols_ =  nvar_ = 2;
    
    init();
}

ViscosityJetMatrix::ViscosityJetMatrix(int nrows, int ncols, int nvar) {
    nrows_ = nrows;
    ncols_ = ncols;
    nvar_  = nvar;
    
    init();
}

ViscosityJetMatrix::~ViscosityJetMatrix(){
}


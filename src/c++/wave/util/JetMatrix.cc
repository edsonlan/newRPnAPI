#include "JetMatrix.h"

JetMatrix::JetMatrix(){
    vec = 0;
}

JetMatrix::JetMatrix(int n){
    vec = 0;

    resize(n);
}

JetMatrix::JetMatrix(int vars, int eqs){
    vec = 0;

    resize(vars, eqs);
}

JetMatrix::JetMatrix(const JetMatrix *o){
    vec = 0;

    copy(*o);
}

JetMatrix::JetMatrix(const JetMatrix &o){
    vec = 0;

    copy(o);
}

JetMatrix::~JetMatrix(){
    if (vec != 0) delete [] vec;
}

std::ostream & operator<<(std::ostream &out, const JetMatrix &r) {
    out << "Function = " << r.function() << "\n\n";
    out << "Jacobian =\n" << r.Jacobian() << "\n\n";

    std::vector<DoubleMatrix> Hessian(r.Hessian());
    for (int i = 0; i < Hessian.size(); i++) out << "Hessian[" << i << "] =\n" << Hessian[i] << "\n\n";

    return out;
}


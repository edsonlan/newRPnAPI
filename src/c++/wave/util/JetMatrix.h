#ifndef _JETMATRIX_
#define _JETMATRIX_

#include <vector>
#include <iostream>
#include "RealVector.h"
#include "DoubleMatrix.h"

class JetMatrix {
    private:
    protected:
        // Number of conservation equations and number of variables.
        //
        int Ceqs, vars;

        // Total number of stored elements (size and space).
        //
        long int vec_size;
        double *vec;

        void resize();
        void copy(const JetMatrix &o);
    public:
        JetMatrix();
        JetMatrix(int n);
        JetMatrix(int v, int e);
        JetMatrix(const JetMatrix *o);
        JetMatrix(const JetMatrix &o);

        ~JetMatrix();

        JetMatrix operator=(const JetMatrix &original);

        void resize(int n);
        void resize(int v, int e);

        int number_of_conservation_equations() const;
        int number_of_variables() const;
        int size() const;

        // Extract the function, Jacobian and Hessian.
        //
        RealVector function() const;
        DoubleMatrix Jacobian() const;
        std::vector<DoubleMatrix> Hessian() const;

        // Extract a matrix from the Hessian, associated to the the i-th component of the mapping's domain.
        // The matrix has components (j, k) of the form:
        //
        //     d^2 F_j / d u_i d u_k
        //
        // where F_j are the functions that form the mapping.
        //
        DoubleMatrix extract_matrix_from_Hessian(int i) const;

        // Set/get the values of the function.
        //
        void set(int i, double x);
        double get(int i) const;

        // Set/get the values of the Jacobian.
        //
        void set(int i, int j, double x);
        double get(int i, int j) const;

        // Set/get the values of the Hessian.
        //
        void set(int i, int j, int k, double x);
        double get(int i, int j, int k) const;

        // Access the data directly
        double * data();
        const double * data() const;

        friend std::ostream & operator<<(std::ostream &out, const JetMatrix &r);
};

inline void JetMatrix::resize(){
    if (vec != 0) delete [] vec;

    vec_size = Ceqs + Ceqs*vars + Ceqs*vars*vars;
    vec = new double[vec_size];
    for (int i = 0; i < vec_size; i++) vec[i] = 0.0;

    return;
}

inline void JetMatrix::copy(const JetMatrix &o){
    Ceqs = o.Ceqs;
    vars = o.vars;

    resize();

    for (int i = 0; i < vec_size; i++) vec[i] = o.vec[i];

    return;
}

inline JetMatrix JetMatrix::operator=(const JetMatrix &original){
    if (this != &original) copy(original);

    return *this;
}

inline void JetMatrix::resize(int n){
    resize(n, n);

    return;
}

inline void JetMatrix::resize(int v, int e){
    vars = v;
    Ceqs = e;

    resize();

    return;
}

inline int JetMatrix::number_of_conservation_equations() const {
    return Ceqs;
}

inline int JetMatrix::number_of_variables() const {
    return vars;
}

inline int JetMatrix::size() const {
    return vec_size;
}

inline RealVector JetMatrix::function() const {
    RealVector f(Ceqs);
    for (int i = 0; i < Ceqs; i++) f(i) = vec[i];

    return f;
}

inline DoubleMatrix JetMatrix::Jacobian() const {
    DoubleMatrix J(Ceqs, vars);
    for (int i = 0; i < Ceqs*vars; i++) J(i) = vec[Ceqs + i];

    return J;
}

inline std::vector<DoubleMatrix> JetMatrix::Hessian() const {
    std::vector<DoubleMatrix> H(Ceqs);
    
    for (int i = 0; i < Ceqs; i++){
        H[i].resize(vars, vars);
        for (int j = 0; j < vars; j++){
            for (int k = 0; k < vars; k++) H[i](j, k) = vec[Ceqs + Ceqs*vars + i*(vars*vars) + j*vars + k];
        }
    }

    return H;
}

// Extract a matrix from the Hessian, associated to the the i-th component of the mapping's domain.
// The matrix has components (j, k) of the form:
//
//     d^2 F_j / d u_i d u_k
//
// where F_j are the functions that form the mapping.
//
inline DoubleMatrix JetMatrix::extract_matrix_from_Hessian(int i) const {
    DoubleMatrix R(Ceqs, vars);

    for (int j = 0; j < Ceqs; j++){
        for (int k = 0; k < vars; k++) R(j, k) = get(j, i, k);
    }

    return R;
}

// Set/get the values of the function.
//
inline void JetMatrix::set(int i, double x){
    vec[i] = x;

    return;
}

inline double JetMatrix::get(int i) const {
    return vec[i];
}

// Set/get the values of the Jacobian.
//
inline void JetMatrix::set(int i, int j, double x){
    vec[Ceqs + i*vars + j] = x;

    return;
}

inline double JetMatrix::get(int i, int j) const {
    return vec[Ceqs + i*vars + j];
}

// Set/get the values of the Hessian.
//
inline void JetMatrix::set(int i, int j, int k, double x){
    vec[Ceqs + Ceqs*vars + i*(vars*vars) + j*vars + k] = x;

    return;
}

inline double JetMatrix::get(int i, int j, int k) const {
    return vec[Ceqs + Ceqs*vars + i*(vars*vars) + j*vars + k];
}

inline double * JetMatrix::data(){
    return vec;
}

inline const double * JetMatrix::data() const {
    return vec;
}

#endif // _JETMATRIX_


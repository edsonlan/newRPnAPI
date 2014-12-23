/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) HessianMatrix.h
 */

#ifndef _HessianMatrix_H
#define _HessianMatrix_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
//#include "Vector.h"
#include "except.h"
#include "RealVector.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

/*! @brief Utility class to store a second derivative matrix
 * 
 *  This matrix has n rows , n columns and deep n
 *
 * @ingroup wave
 */


class HessianMatrix {
private:
    int n_comps_, size_;
    RealVector v_;

    class RangeViolation : public exception {
    };
public:
    HessianMatrix(void);
    HessianMatrix(const int n_comps);
    HessianMatrix(const HessianMatrix & hessianMatrix);

    int n_comps(void) const;
    void resize(int n_comps);
    void range_check(int comp) const;
    HessianMatrix &zero(void);

    double * operator()(void);
    double operator()(int i, int j, int k) const;
    void operator()(int i, int j, int k, double value);

};

inline HessianMatrix::HessianMatrix(void) :
n_comps_(1),
size_(1),
v_(1) {
}

inline HessianMatrix::HessianMatrix(const int n_comps) :
n_comps_(n_comps),
size_(n_comps * n_comps * n_comps),
v_(size_) {
}

inline HessianMatrix::HessianMatrix(const HessianMatrix & hessianMatrix) :
n_comps_(hessianMatrix.n_comps_),
size_(n_comps_ * n_comps_ * n_comps_),
v_(hessianMatrix.v_) {
}

inline int HessianMatrix::n_comps(void) const {
    return n_comps_;
}

inline void HessianMatrix::resize(int n_comps) {
    v_.resize(n_comps * n_comps * n_comps);
    n_comps_ = n_comps;
}

inline void HessianMatrix::range_check(int comp) const {
    if (comp < 0 || comp >= n_comps())
        throw (HessianMatrix::RangeViolation());
}

inline HessianMatrix & HessianMatrix::zero(void) {
    for (int i = 0; i < v_.size(); i++) {
        v_(i) = 0;
    }


    return *this;
}

inline double * HessianMatrix::operator()(void) {
    return v_.components();
}

inline double HessianMatrix::operator()(int i, int j, int k) const {
    range_check(i);
    range_check(j);
    range_check(k);
    return v_.component(i * n_comps_ * n_comps_ + j * n_comps_ + k);
}

inline void HessianMatrix::operator()(int i, int j, int k, double value) {
    range_check(i);
    range_check(j);
    range_check(k);
    double * value_ = &v_.component(i * n_comps_ * n_comps_ + j * n_comps_ + k);
    *value_ = value;
}



#endif //! _HessianMatrix_H

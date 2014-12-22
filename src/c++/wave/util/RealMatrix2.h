/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RealMatrix2.h
 */

#ifndef _RealMatrix2_H
#define _RealMatrix2_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */

#include "RealVector.h"
#include "eigen.h"

//#include "Vector.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


extern"C" {
    void dgetrf_(int *, int *, double *, int *, int *, int *);
}

class RealMatrix2 {
private:

    RealVector * data_;
    int row_, col_;

      class RangeViolation : public exception {
    };

public:

    RealMatrix2(int, int);
    RealMatrix2(int);

    virtual ~RealMatrix2();
    RealMatrix2();
    RealMatrix2(const RealMatrix2 &);

    void range_check(int i, int j) const;

    void resize(int row, int col);

    double operator ()(int i, int j) const;

    void operator ()(int i, int j, double value);

    RealMatrix2 & operator=(const RealMatrix2 &);
    
    bool operator==(const RealMatrix2 &);

    bool operator != (const RealMatrix2 &);

    RealMatrix2 & operator-(const RealMatrix2 &);

    RealMatrix2 & operator+(const RealMatrix2 &);

    void getColumn(int, double *);

    void setColumn(int, double *);

    void copySubMatrix(const int, const int, const int, const int, const int, const int, RealMatrix2 &);

    double determinant();

    void mul(const RealMatrix2 &)const;

    void fillEigenData(int stateSpaceDim, RealMatrix2 & df, double & eigenValR, double & eigenValI, RealVector & eigenVec);

    void scale(double t);

    friend std::ostream & operator<<(std::ostream&, const RealMatrix2&);

    RealMatrix2 & zero();

    int row()const;

    int col()const;

};

inline void RealMatrix2::resize(int row, int col) {
    row_ = row;
    col_ = col;
    data_->resize(row * col);
}

inline void RealMatrix2::range_check(int i, int j) const {
    if (((i < 0) && (i >= row_)) || ((j < 0) && (j >= col_)))
        throw(RealMatrix2::RangeViolation());
}

inline double RealMatrix2::operator()(int i, int j) const {
    range_check(i, j);
    return data_->component(i * row_ + j);

}

inline RealMatrix2 & RealMatrix2::operator=(const RealMatrix2 & source) {
    
    if ((source.row_!=row_) || (source.col_!=col_))
         throw(RealMatrix2::RangeViolation());
    
    for (int i=0;i < row_;i++){
        for (int j=0;j< col_;j++){
            operator()(i,j,source(i,j));
        }
    }
    
    return *this;
    
    
}


inline RealMatrix2 & RealMatrix2::operator-(const RealMatrix2 & b) {
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            operator()(i, j, operator()(i, j) - b(i, j));
        }
    }
    return *this;

}

inline bool RealMatrix2::operator==(const RealMatrix2& matrix) {

    if ((matrix.row() != row_) || (matrix.col() != col_))
        return false;

    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            if (matrix(i, j) != operator()(i, j))
                return false;
        }
    }

    return true;

}


inline bool RealMatrix2::operator!=(const RealMatrix2& matrix) {

    
    return (!(*this==matrix));
    
//    if ((matrix.row() != row_) || (matrix.col() != col_))
//        return true;
//
//    for (int i = 0; i < row_; i++) {
//        for (int j = 0; j < col_; j++) {
//            if (matrix(i, j) != operator()(i, j))
//                return true;
//        }
//    }
//
//    return false;

}




inline RealMatrix2 & RealMatrix2::operator+(const RealMatrix2 & b) {
    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            operator()(i, j, operator()(i, j) + b(i, j));
        }
    }
    return *this;

}

inline RealMatrix2 & RealMatrix2::zero(void) {
    for (int i = 0; i < data_->size(); i++) {
        data_->component(i)=0;
    }
    return *this;
}

inline void RealMatrix2::operator ()(int i, int j, double value) {
    range_check(i, j);
    data_->component(i * row_ + j) = value;

}

inline int RealMatrix2::row()const {
    return row_;
}

inline int RealMatrix2::col()const {
    return col_;
}

inline void RealMatrix2::scale(double t) {

    for (int i = 0; i < row_; i++) {
        for (int j = 0; j < col_; j++) {
            data_->component(i * row_ + j) = data_->component(i * row_ + j) * t;
        }
    }

}
#endif //! _RealMatrix2_H

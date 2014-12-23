#ifndef _DOUBLEMATRIX_
#define _DOUBLEMATRIX_

#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <limits>
#include "Matrix.h"

#define DOUBLEMATRIXPRINTWIDTH 14

extern "C" {
    // Linear system of equations solver
    void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);

    // LU decomposition
    void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

    // Matrix inversion
    void dgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);

    // QR decomposition
    void dgeqp3_(int *M, int *N, double *A, int *LDA, int *JPVT, double *TAU, double *WORK, int *LWORK, int *INFO);
}

class DoubleMatrix : public Matrix<double> {
    private:
    protected:
        int w_; // Default width for the centered_text.
        static std::string centered_text(double x, int w);

        // For Gram-Schmidt:
        //
        double inner_product_by_rows(int row1, int row2);
        void proj_u_v_by_rows(int u_row, int v_row, double *w);

        double inner_product_by_columns(int col1, int col2);
        void proj_u_v_by_columns(int u_col, int v_col, double *w);
    public:
        DoubleMatrix(void);
        DoubleMatrix(int n, int m);
        DoubleMatrix(const DoubleMatrix &original);
        DoubleMatrix(const DoubleMatrix *original);
        DoubleMatrix(const Matrix<double> &original);
        DoubleMatrix(const Matrix<double> *original);
        DoubleMatrix(int n, int m, const double *original);

        virtual ~DoubleMatrix();

        DoubleMatrix operator=(const DoubleMatrix &original);

        // Deprecated:
        // void print(void) const;

        // Create the identity matrix.
        //
        static DoubleMatrix eye(int n);

        // Create a zero-filled matrix.
        //
        static DoubleMatrix zero(int m, int n);
        
        // Maximum and minimum.
        //
        double max();
        double min();
        void minmax(double &minimum, double &maximum);

        double norm_max();

        // Set/get the width used for printing
        void w(int ww){w_ = ww; return;}
        int w(void) const {return w_;}

        // Gram-Schmidt orthogonalization process.
        void Gram_Schmidt_orthogonalization_by_rows();
        void Gram_Schmidt_orthogonalization_by_columns();

        // Extraction of rows and columns
        DoubleMatrix row(int r);
        DoubleMatrix column(int c);
        
        // Matrix sum and subtraction
        friend DoubleMatrix sum(const DoubleMatrix &A, const DoubleMatrix &B);
        friend DoubleMatrix operator+(const DoubleMatrix &A, const DoubleMatrix &B);

        friend DoubleMatrix sub(const DoubleMatrix &A, const DoubleMatrix &B);
        friend DoubleMatrix operator-(const DoubleMatrix &A, const DoubleMatrix &B);

        // Matrix multiplication 
        friend DoubleMatrix mult(const DoubleMatrix &A, const DoubleMatrix &B);
        friend DoubleMatrix operator*(const DoubleMatrix &A, const DoubleMatrix &B);

        // Multiplication by a scalar
        friend DoubleMatrix mult(const DoubleMatrix &A, double alpha);
        friend DoubleMatrix mult(double alpha, const DoubleMatrix &A);
        friend DoubleMatrix operator*(const DoubleMatrix &A, double alpha);
        friend DoubleMatrix operator*(double alpha, const DoubleMatrix &A);

        // Some matrix operations
        friend int solve(const DoubleMatrix &A, const DoubleMatrix &b, DoubleMatrix &x);
        friend int inverse(const DoubleMatrix &A, DoubleMatrix &B);
        friend DoubleMatrix inverse(const DoubleMatrix &A);

        friend double det(const DoubleMatrix &A);
        friend double derivative_det(const DoubleMatrix &A, const DoubleMatrix &derivative_A);

        friend DoubleMatrix transpose(const DoubleMatrix &A);
        friend int QR(const DoubleMatrix &A, DoubleMatrix &Q, DoubleMatrix &R);
  
        // Trace
        friend double trace(const DoubleMatrix &A);

        // Output to a stream
        friend std::ostream& operator<<(std::ostream& stream, const DoubleMatrix &m); 
};

#endif // _DOUBLEMATRIX_


#ifndef _MATRIX_
#define _MATRIX_

#include <stdlib.h>
#include <cmath>
#include <iostream>

/*    SIMPLE USE EXAMPLE:

#include "Matrix.h"
#include <stdio.h>

int main(){
    int n = 2, m = 3;
    Matrix<double> A(n, m);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) A(i, j) = 1 + 2*i + 3*j;
    }

    Matrix<double> B(A);

    Matrix<double> C, D;
    C = D = B;

    printf("C.rows = %d, C.cols = %d\n", C.rows(), C.cols());

    for (int i = 0; i < C.rows(); i++){
        for (int j = 0; j < C.cols(); j++) printf("C(%d, %d) = %f\n", i, j, C(i, j));
    }

    return 0;
}

     END OF EXAMPLE */

template <typename T>
class Matrix {
    private:

    protected:
        T *vec;
        int rows_, cols_;

        int min(int x, int y);

        void copy(int n, int m, const T *orig);

    public:
        Matrix();
        Matrix(int n, int m);
        Matrix(const Matrix<T> &original);
        Matrix(const Matrix<T> *original);
        Matrix(int n, int m, const T *original);
        virtual ~Matrix();

        virtual T&       operator()(int i, int j);
        virtual T const& operator()(int i, int j) const;

        virtual T&       operator()(int i);
        virtual T const& operator()(int i) const;

        virtual T*       data(void);
        virtual T const* data(void) const;

        Matrix<T> operator=(const Matrix<T> &original);

        void resize(int newn, int newm);

        int rows(void) const;
        int cols(void) const;

        void size(int &r, int &c) const;
        
        // Extract submatrix.
        //
        Matrix<T> submatrix(int minr, int maxr, int minc, int maxc);
        Matrix<T> submatrix(int minr, int maxr, int minc, int maxc) const;

        friend std::ostream& operator<<(std::ostream& stream, const Matrix<T> &m){
            for (int i = 0; i < m.rows_; i++){
                stream << "|";
                for (int j = 0; j < m.cols_; j++){
                    //stream << DoubleMatrix::centered_text(m(i, j), m.w_);
                    stream << " " << m(i, j) << " ";
                    if (j < m.cols_ - 1) stream << ":";
                }
                stream << "|" << std::endl;
            }

            return stream;
        }
};

template <typename T> int Matrix<T>::min(int x, int y){
    return (x < y) ? x : y;
}

template <typename T> void Matrix<T>::copy(int n, int m, const T *orig){
    resize(n, m);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n*m; i++) vec[i] = orig[i];

    return;
}

template <typename T> Matrix<T>::Matrix(){
    vec = 0;

    cols_ = rows_ = 0;
    
    //std::cout << "Matrix<T> 1, this = " << this << std::endl;
}

template <typename T> Matrix<T>::Matrix(int n, int m){
    vec = 0;

    resize(n, m);
    
    //std::cout << "Matrix<T> 2, this = " << this << std::endl;
}

template <typename T> Matrix<T>::Matrix(const Matrix<T> &original){
    vec = 0;

    copy(original.rows_, original.cols_, original.vec);
    
    //std::cout << "Matrix<T> 3, this = " << this << std::endl;
}

template <typename T> Matrix<T>::Matrix(const Matrix<T> *original){
    vec = 0;

    copy(original->rows_, original->cols_, original->vec);
    
    //std::cout << "Matrix<T> 4, this = " << this << std::endl;
}

template <typename T> Matrix<T>::Matrix(int n, int m, const T *original){
    vec = 0;

    copy(n, m, original);
    
    //std::cout << "Matrix<T> 5, this = " << this << std::endl;
}

template <typename T> Matrix<T>::~Matrix(){
    if (vec != 0) delete [] vec;
    //std::cout << "Matrix<T> dtor, this = " << this << std::endl;
}

template <typename T> T& Matrix<T>::operator()(int i, int j){
    return vec[i*cols_ + j];
}

template <typename T> T const& Matrix<T>::operator()(int i, int j) const {
    return vec[i*cols_ + j];
}

template <typename T> T& Matrix<T>::operator()(int i){
    return vec[i];
}

template <typename T> T const& Matrix<T>::operator()(int i) const {
    return vec[i];
}

template <typename T> T* Matrix<T>::data(void){
    return &vec[0];
}

template <typename T> T const* Matrix<T>::data(void) const {
    return &vec[0];
}

template <typename T> Matrix<T> Matrix<T>::operator=(const Matrix<T> &original){
    if (this != &original) copy(original.rows_, original.cols_, original.vec);

    return *this;
}

template <typename T> void Matrix<T>::resize(int newn, int newm){
    if (newn == 0 || newm == 0) return;

    if (vec == 0) vec = new T[newn*newm];
    else {
        T *temp = new T[newn * newm];

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < std::min(rows_, newn); i++){
            for (int j = 0; j < std::min(cols_, newm); j++){
                temp[i*newm + j] = vec[i*cols_ + j];
            }
        }

        delete [] vec;

        vec = temp;
    }

    rows_ = newn;
    cols_ = newm;

    return;
}

template <typename T> int Matrix<T>::rows(void) const {
    return rows_;
}

template <typename T> int Matrix<T>::cols(void) const {
    return cols_;
}

template <typename T> void Matrix<T>::size(int &r, int &c) const {
    r = rows_; 
    c = cols_; 
    return;
}

template <typename T> Matrix<T> Matrix<T>::submatrix(int minr, int maxr, int minc, int maxc){
	Matrix<T> subM(maxr - minr + 1, maxc - minc + 1);
	
	for (int i = 0; i < subM.rows(); i++){
		for (int j = 0; j < subM.cols(); j++){
		    subM(i, j) = this->operator()(i + minr, j + minc);
		}
	}
	
	return subM;
}

template <typename T> Matrix<T> Matrix<T>::submatrix(int minr, int maxr, int minc, int maxc) const {
	Matrix<T> subM(maxr - minr + 1, maxc - minc + 1);
	
	for (int i = 0; i < subM.rows(); i++){
		for (int j = 0; j < subM.cols(); j++){
		    subM(i, j) = this->operator()(i + minr, j + minc);
		}
	}
	
	return subM;
}

#endif // _MATRIX_


/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RealMatrix2.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "RealMatrix2.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */





RealMatrix2::RealMatrix2(int row, int col) : data_(new RealVector(row*col)), row_(row), col_(col) {
}

RealMatrix2::RealMatrix2(int rowcol) : data_(new RealVector(rowcol*rowcol)), row_(rowcol), col_(rowcol) {

    for (int i = 0; i < rowcol; i++) {
        for (int j = 0; j < rowcol; j++) {


            if (i == j)
                data_->component(i * rowcol + j) = 1;
            else
                data_->component(i * rowcol + j) = 0;
        }
    }

}

RealMatrix2::RealMatrix2() : data_(new RealVector(4)), row_(2), col_(2) {
}

RealMatrix2::RealMatrix2(const RealMatrix2 & copy) : data_(copy.data_) {
}

void RealMatrix2::fillEigenData(int stateSpaceDim, RealMatrix2 &df, double & eigenValR, double & eigenValI, RealVector & eigenVec) {


    int i, j, info, lwork = 5 * stateSpaceDim;

    int lda = stateSpaceDim;
    int ldvr = stateSpaceDim;
    int ldvl = stateSpaceDim;

    double wr[stateSpaceDim];
    double wi[stateSpaceDim];

    double J[stateSpaceDim][stateSpaceDim], vl[stateSpaceDim][stateSpaceDim], vr[stateSpaceDim][stateSpaceDim];

    double work[stateSpaceDim];

    for (i = 0; i < stateSpaceDim; i++) {
        for (j = 0; j < stateSpaceDim; j++) {
            J[i][j] = df(i, j);
        }
    }

    dgeev_("N", "V", &stateSpaceDim, &J[0][0], &lda, &wr[0], &wi[0],
            &vl[0][0], &ldvl, &vr[0][0], &ldvr, &work[0], &lwork,
            &info); // only the right eigenvectors are needed.

}

double RealMatrix2::determinant() {

    int numRow = row();
    int numCol = col();
    int INFO;
    double result = 1;

    if ((numRow == 2) && (numCol == 2))
        result = operator()(0, 0) * operator()(1, 1) - operator()(0, 1) * operator()(1, 0);
    else
        if ((numRow == 3) && (numCol == 3))
        result = operator()(0, 0) * operator()(1, 1) * operator()(2, 2)
        + operator()(0, 1) * operator()(1, 2) * operator()(2, 0)
        + operator()(0, 2) * operator()(1, 0) * operator()(2, 1)
        - operator()(0, 2) * operator()(1, 1) * operator()(2, 0)
        - operator()(0, 1) * operator()(1, 0) * operator()(2, 2)
        - operator()(0, 0) * operator()(1, 2) * operator()(2, 1);
    else {


        double A[numRow][numCol];
        for (int i = 0; i < numRow; i++) {
            for (int j = 0; j < numCol; j++) {
                A[i][j] = operator()(i, j);
            }
        }

        int * IPIV;
        int IPIV_length;

        if (numRow < numCol) {
            IPIV = new int[numRow];
            IPIV_length = numRow;
        } else {
            IPIV = new int[numCol];
            IPIV_length = numCol;

        }

        dgetrf_(&numRow, &numCol, &A[0][0], &numRow, IPIV, &INFO);

        RealMatrix2 LU(numRow, numCol);

        for (int i = 0; i < numRow; i++) {
            for (int j = 0; j < numCol; j++) {
                LU(i, j, A[i][j]);
            }
        }

        for (int i = 0; i < LU.col(); i++)
            result *= LU(i, i);

        for (int i = 0; i < IPIV_length; i++)
            if (IPIV[i] != i + 1)
                result = -result;
        delete [] IPIV;
    }

    return result;
}

void RealMatrix2::mul(const RealMatrix2 & m1) const {

    int i, j, k;

    if (col_ != m1.row() || col_ != m1.col())

        throw(RealMatrix2::RangeViolation());

    RealVector temp(row_ * col_);

    for (i = 0; i < row_; i++) {
        for (j = 0; j < col_; j++) {
            temp.component((row_) + (i * col_ + j)) = 0.0;
            for (k = 0; k < col_; k++) {
                temp.component((row_) + (i * col_ + j)) += data_->component((row_) + (i * col_ + k)) * m1(k, j);
            }
        }
    }

    for (i = 0; i < row_; i++) {
        for (j = 0; j < col_; j++) {
            data_->component((row_) + (i * col_ + j)) = temp((row_) + (i * col_ + j));
        }
    }
}

void RealMatrix2::getColumn(int colIndex, double * outputArray) {
    if (colIndex < 0 || colIndex >= col_)
        throw(RealMatrix2::RangeViolation());

    for (int i = 0; i < row_; i++) {
        outputArray[i] = operator()(i, colIndex);
    }

}

void RealMatrix2::setColumn(int colIndex, double * inputArray) {
    if (colIndex < 0 || colIndex >= col_)
        throw(RealMatrix2::RangeViolation());

    for (int i = 0; i < row_; i++) {
        operator()(i, colIndex, inputArray[i]);
    }
}

void RealMatrix2::copySubMatrix(const int rowSource, const int colSource,
        const int numRow, const int numCol, const int rowDest,
        const int colDest, RealMatrix2 & target) {
    int i, j;

    if (*this != target) {
        for (i = 0; i < numRow; i++) {
            for (j = 0; j < numCol; j++) {
                target(rowDest + i, colDest + j, operator()(rowSource + i, colSource + j));
            }
        }
    } else {
        double tmp[numRow][numCol];
        for (i = 0; i < numRow; i++) {
            for (j = 0; j < numCol; j++) {
                tmp[i][j] = operator()(rowSource + i, colSource + j);
            }
        }
        for (i = 0; i < numRow; i++) {
            for (j = 0; j < numCol; j++) {
                target(rowDest + i, colDest + j, tmp[i][j]);
            }
        }
    }
}

std::ostream & operator<<(std::ostream& o, const RealMatrix2& jMatrix) {

    for (int i = 0; i < jMatrix.row_; i++) {
        for (int j = 0; j < jMatrix.col_; j++) {

            o << "(" << i << "," << j << ")" << jMatrix(i, j) << "\n";
        }
    }

    return o;
}

RealMatrix2::~RealMatrix2() {
    delete data_;
}


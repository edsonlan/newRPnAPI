#include "DoubleMatrix.h"

std::string DoubleMatrix::centered_text(double x, int w){
    std::stringstream s;
    s << x;

    std::string temp(s.str());

    int noel, sl, noer;
    sl = temp.size();
    noel = floor((w - sl)/2);
    noer = w - sl - noel;

    std::stringstream ss;
    for (int i = 0; i < noel; i++) ss << std::string(" ");

    ss << x;

    for (int i = 0; i < noer; i++) ss << std::string(" ");

    return ss.str();
}

// For the Gram-Schmidt method: inner product of two rows.
//
double DoubleMatrix::inner_product_by_rows(int row1, int row2){
    double p = 0.0;

    for (int j = 0; j < cols_; j++) p += vec[row1*cols_ + j]*vec[row2*cols_ + j];

    return p;
}

// For the Gram-Schmidt method: projection of one row over another one.
//
void DoubleMatrix::proj_u_v_by_rows(int u_row, int v_row, double *w){
    double p = inner_product_by_rows(v_row, u_row)/inner_product_by_rows(u_row, u_row);

    for (int i = 0; i < cols_; i++) w[i] = p*vec[u_row*cols_ + i];

    return;
}

// For the Gram-Schmidt method: inner product of two cols.
//
double DoubleMatrix::inner_product_by_columns(int col1, int col2){
    double p = 0.0;

    for (int j = 0; j < rows_; j++) p += vec[j*cols_ + col1]*vec[j*cols_ + col2];

    return p;
}

// For the Gram-Schmidt method: projection of one column over another one.
//
void DoubleMatrix::proj_u_v_by_columns(int u_col, int v_col, double *w){
    double p = inner_product_by_columns(v_col, u_col)/inner_product_by_columns(u_col, u_col);

    for (int i = 0; i < rows_; i++) w[i] = p*vec[i*cols_ + u_col];

    return;
}

// Extraction of rows and columns
DoubleMatrix DoubleMatrix::row(int r){
    DoubleMatrix rr(1, cols_);

    for (int i = 0; i < cols_; i++) rr(0, i) = vec[r*cols_ + i];

    return rr;
}

DoubleMatrix DoubleMatrix::column(int c){
    DoubleMatrix cc(rows_, 1);

    for (int i = 0; i < rows_; i++) cc(i, 0) = vec[i*cols_ + c];

    return cc;
}

DoubleMatrix::DoubleMatrix(void) : Matrix<double>(), w_(DOUBLEMATRIXPRINTWIDTH){
    //std::cout << "DM 1, this = " << this << std::endl;
}

DoubleMatrix::DoubleMatrix(int n, int m) : Matrix<double>(n, m), w_(DOUBLEMATRIXPRINTWIDTH){
    //std::cout << "DM 2, this = " << this << std::endl;
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix &original) : Matrix<double>(original), w_(DOUBLEMATRIXPRINTWIDTH){
    //std::cout << "DM 3, this = " << this << std::endl;
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix *original) : Matrix<double>(original), w_(DOUBLEMATRIXPRINTWIDTH){
    //std::cout << "DM 4, this = " << this << std::endl;
}

DoubleMatrix::DoubleMatrix(const Matrix<double> &original) : Matrix<double>(original), w_(DOUBLEMATRIXPRINTWIDTH){
    //std::cout << "DM 5, this = " << this << std::endl;
}

DoubleMatrix::DoubleMatrix(const Matrix<double> *original) : Matrix<double>(original), w_(DOUBLEMATRIXPRINTWIDTH){
    //std::cout << "DM 6, this = " << this << std::endl;
}

DoubleMatrix::DoubleMatrix(int n, int m, const double *original) : Matrix<double>(n, m, original), w_(DOUBLEMATRIXPRINTWIDTH){
    //std::cout << "DM 7, this = " << this << std::endl;
}

DoubleMatrix::~DoubleMatrix(){
    //std::cout << "DM dtor, this = " << this << std::endl;
}

// Deprecated.
//
//void DoubleMatrix::print(void) const {
//    for (int i = 0; i < rows_; i++){
//        printf("|");
//        for (int j = 0; j < cols_; j++){
//            //printf(" % 14.8g ", vec[i*cols_ + j]);
//            printf("%s", centered_text(vec[i*cols_ + j], w_).c_str());
//            if (j < cols_ - 1) printf(":");
//        }
//        printf("|\n");
//    }

//    return;
//}

DoubleMatrix DoubleMatrix::operator=(const DoubleMatrix &original){
    if (this != &original) copy(original.rows_, original.cols_, original.vec);

    return *this;
}

// Identity matrix
//
DoubleMatrix DoubleMatrix::eye(int n){
    DoubleMatrix e(n, n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++){
        for (int j = 0; j < i; j++) e(i, j) = e(j, i) = 0.0;
        e(i, i) = 1.0;
    }

    return e;
}

DoubleMatrix DoubleMatrix::zero(int m, int n){
    DoubleMatrix z(m, n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++) z(i, j) = 0.0;
    }

    return z;
}


double DoubleMatrix::max(){
	double minimum, maximum;
	
	minmax(minimum, maximum);
	
	return maximum;
}

double DoubleMatrix::min(){
	double minimum, maximum;
	
	minmax(minimum, maximum);
	
	return minimum;
}

void DoubleMatrix::minmax(double &minimum, double &maximum){
    minimum =  std::numeric_limits<double>::max();
    maximum = -std::numeric_limits<double>::max();
	
    for (int i = 0; i < rows_*cols_; i++){
        if (vec[i] > maximum) maximum = vec[i];
        if (vec[i] < minimum) minimum = vec[i];
    }
	
    return;
}

double DoubleMatrix::norm_max(){
    double d = 0.0;

    for (int i = 0; i < rows_*cols_; i++) d = std::max(d, std::abs(vec[i]));

    return d;
}

// Returns A + B
//
DoubleMatrix sum(const DoubleMatrix &A, const DoubleMatrix &B){
    DoubleMatrix C(A.rows_, A.cols_);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < A.rows_*A.cols_; i++) C(i) = A(i) + B(i);

    return C;
}

DoubleMatrix operator+(const DoubleMatrix &A, const DoubleMatrix &B){
    return sum(A, B);
}

// Returns A - B
//
DoubleMatrix sub(const DoubleMatrix &A, const DoubleMatrix &B){
    DoubleMatrix C(A.rows_, A.cols_);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < A.rows_*A.cols_; i++) C(i) = A(i) - B(i);

    return C;
}

DoubleMatrix operator-(const DoubleMatrix &A, const DoubleMatrix &B){
    return sub(A, B);
}

// Returns A*B
//
DoubleMatrix mult(const DoubleMatrix &A, const DoubleMatrix &B){
    int m = A.rows_;
    int n = B.cols_;
    int p = A.cols_;

    DoubleMatrix C(m, n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            double s = 0.0;

            for (int k = 0; k < p; k++) s += A(i, k)*B(k, j);
            C(i, j) = s;
        }
    }

    return C;
}

DoubleMatrix operator*(const DoubleMatrix &A, const DoubleMatrix &B){
    return mult(A, B);
}

// Returns alpha*A
//
DoubleMatrix mult(const DoubleMatrix &A, double alpha){
    DoubleMatrix B(A);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < B.rows_*B.cols_; i++) B(i) *= alpha;

    return B;
}

DoubleMatrix mult(double alpha, const DoubleMatrix &A){
    return mult(A, alpha);
}

DoubleMatrix operator*(const DoubleMatrix &A, double alpha){
    return mult(A, alpha);
}

DoubleMatrix operator*(double alpha, const DoubleMatrix &A){
    return mult(A, alpha);
}

int solve(const DoubleMatrix &A, const DoubleMatrix &b, DoubleMatrix &x){
    int i, j;
    int dim = A.rows();
    int nrhs = 1;
    int lda = dim;
    int ipiv[dim];
    int ldb = dim;
    int info;

    // Create a transposed copy of A to be used by LAPACK's dgesv:
    DoubleMatrix B(transpose(A));

    // Create a copy of b to be used by LAPACK's dgesv:
    double bb[dim];
    for (i = 0; i < dim; i++) bb[i] = b(i, 0);

    dgesv_(&dim, &nrhs, B.data(), &lda, &ipiv[0], bb, &ldb, &info);

    x.resize(dim, 1);
    if (info == 0){
        for (i = 0; i < dim; i++) x(i, 0) = bb[i];
    }

    return info;
}

int inverse(const DoubleMatrix &A, DoubleMatrix &B){
    int n = A.rows_;

    B = transpose(A);

    // LU factorization
    //
    int lu_info, lda = n, ipiv[n];
    dgetrf_(&n, &n, B.data(), &lda, ipiv, &lu_info);
    if (lu_info != 0) return lu_info;

    //for (int i = 0; i < n; i++) printf("ipiv(%d) = %d\n", i, ipiv[i]);

    // Matrix inversion proper
    //
    int lwork = n;
    double work[lwork];
    int inv_info;
    dgetri_(&n, B.data(), &lda, ipiv, work, &lwork, &inv_info);

    B = transpose(B);

    return inv_info;
}

// Only to be used if matrix A is known to be non-singular.
//
DoubleMatrix inverse(const DoubleMatrix &A){
    DoubleMatrix B; 
    inverse(A, B); 

    return B;
}

//double det(const DoubleMatrix &A){
//    int n = A.rows(); 

//    if (n == 1) return A(0, 0);
//    else {
//        double *det_AA = new double[n];

//        #pragma omp parallel for schedule(dynamic)
//        for (int i = 0; i < n; i++){
//            DoubleMatrix AA(n - 1, n - 1);
//            if (i == 0){
//                for (int j = 1; j < n; j++){
//                    for (int k = 1; k < n; k++) AA(j - 1, k - 1) = A(j, k);
//                }
//            }
//            else if (i == (n - 1)){
//                for (int j = 1; j < n; j++){
//                    for (int k = 0; k < n - 1; k++) AA(j - 1, k - 0) = A(j, k);
//                }
//            }
//            else {
//                for (int j = 1; j < n; j++){
//                    for (int k = 0; k <= i - 1; k++) AA(j - 1, k - 0) = A(j, k);
//                    for (int k = i + 1; k < n; k++) AA(j - 1, k - 1) = A(j, k);
//                }
//            }

//            det_AA[i] = det(AA);
//        }

//        double temp = 0.0, power = 1.0;
//        for (int i = 0; i < n; i++){
//            temp += A(0, i)*det_AA[i]*power;
//            power *= -1.0;
//        }

//        delete [] det_AA;

//        return temp;
//    }
//}

double det(const DoubleMatrix &A){
    int n = A.rows_;

    if (n == 1){
        return A(0, 0);
    }
    else if (n == 2){
        return A(0, 0)*A(1, 1) - A(0, 1)*A(1, 0);
    }
    else if (n == 3) {
        return A(0, 0)*(A(1, 1)*A(2, 2) - A(2, 1)*A(1, 2)) - 
               A(0, 1)*(A(1, 0)*A(2, 2) - A(2, 0)*A(1, 2)) + 
               A(0, 2)*(A(1, 0)*A(2, 1) - A(2, 0)*A(1, 1));
    }
    else {
        // TODO: THIS IS WRONG. Use a handrolled LU for this one.

        DoubleMatrix B(transpose(A));

        // LU factorization
        //
        int lu_info, lda = n, ipiv[n];
        dgetrf_(&n, &n, B.data(), &lda, ipiv, &lu_info);

        if (lu_info != 0) return 0.0;

        double d = 1.0;
        for (int i = 0; i < n; i++){
            if (ipiv[i] == i) d *=  B(i, i);
            else              d *= -B(i, i);
        }

        return d;
    }
}

double derivative_det(const DoubleMatrix &A, const DoubleMatrix &derivative_A){
    double dd = 0.0;

    for (int i = 0; i < A.rows(); i++){
        DoubleMatrix temp(A);

        for (int j = 0; j < A.cols(); j++) temp(i, j) = derivative_A(i, j);

        dd += det(temp);
    }

    return dd;
}

DoubleMatrix transpose(const DoubleMatrix &A){
    DoubleMatrix B(A.cols_, A.rows_);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < A.rows_; i++){
        for (int j = 0; j < A.cols_; j++) B(j, i) = A(i, j);
    }

    return B;
}

// QR decomposition
//
// http://www.netlib.org/lapack/double/dgeqp3.f
//
int QR(const DoubleMatrix &A, DoubleMatrix &Q, DoubleMatrix &R){
    int M = A.rows();
    int N = A.cols();
    int LDA = std::max(1, M);
    int LWORK = 3*N + 1;
    int INFO;
    
    int *JPVT = new int[N]; for (int i = 0; i < N; i++) JPVT[i] = 0.0;
    double *TAU = new double[std::min(M, N)];
    double *WORK = new double[std::max(1, LWORK)];

    DoubleMatrix B = transpose(A); std::cout << "Here. A = \n" << A << std::endl; 

    dgeqp3_(&M, &N, B.data(), &LDA, JPVT, TAU, WORK, &LWORK, &INFO);

    // Output (B is transposed!)
    //
    R.resize(std::max(M, N), N); std::cout << "Here. B = \n" << transpose(B) << std::endl; 

    for (int i = 0; i < R.rows_; i++){
        for (int j = 0; j < i; j++)       R(i, j) = 0.0;
        for (int j = i; j < R.cols_; j++) R(i, j) = B(j, i);
    }

    //Q.resize(M, N);
    Q = DoubleMatrix::eye(M);
    for (int i = 0; i < std::min(M, N); i++){
        DoubleMatrix v(M, 1);
        for (int j = 0; j < i; j++) v(j) = 0.0;
                                    v(i) = 1.0;
        for (int j = i + 1; j < M; j++) v(j) = B(j, i);

        Q = Q*(DoubleMatrix::eye(M) - TAU[i]*v*transpose(v));
    }

    std::cout << "JPVT = (";
    for (int i = 0; i < N; i++) std::cout << " " << JPVT[i];
    std::cout << ")" << std::endl; 

    delete [] WORK;
    delete [] TAU;
    delete [] JPVT;

    return INFO;
}

double trace(const DoubleMatrix &A){
    double t = 0.0;

    for (int i = 0; i < A.rows_; i++) t += A(i, i);

    return t;
}

// Output to a stream
std::ostream& operator<<(std::ostream& stream, const DoubleMatrix &m){
    std::cout << std::endl;
    for (int i = 0; i < m.rows_; i++){
        stream << "|";
        for (int j = 0; j < m.cols_; j++){
            stream << DoubleMatrix::centered_text(m(i, j), m.w_);
            if (j < m.cols_ - 1) stream << ":";
        }
        stream << "|" << std::endl;
    }

    return stream;
}

// Gram-Schmidt orthogonalization process (by rows).
//
void DoubleMatrix::Gram_Schmidt_orthogonalization_by_rows(){
    for (int j = 0; j < rows_; j++){
        for (int i = 0; i < j; i++){
            double temp[cols_];
            proj_u_v_by_rows(i, j, temp);

            for (int k = 0; k < cols_; k++) vec[j*cols_ + k] -= temp[k];
        }
//        normalize(n, &A[j*n]);
    }

    return;
}

// Gram-Schmidt orthogonalization process (by columns).
//
void DoubleMatrix::Gram_Schmidt_orthogonalization_by_columns(){
    double *temp = new double[rows_];

    for (int j = 0; j < cols_; j++){
        for (int i = 0; i < j; i++){

            proj_u_v_by_columns(i, j, temp);

            for (int k = 0; k < rows_; k++) vec[k*cols_ + j] -= temp[k];
        }
//        normalize(n, &A[j*n]);
    }

    delete [] temp;

    return;
}


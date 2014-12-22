#include "RealVector.h"

inline int RealVector::min(int x, int y) {
    return (x < y) ? x : y;
}

RealVector::RealVector(void) {
    data = NULL;
    size_ = 0;
}

RealVector::RealVector(int n) : size_(n) {
    data = new double[n];
}

RealVector::RealVector(int n, const double * values) : size_(n) {
    data = new double[n];
    for (int i = 0; i < n; i++) {
        data[i] = values[i];
    }
}

RealVector::RealVector(int init, int size, double *values) : size_(size){
    data = new double[size];
    
    for (int i = 0; i < size; i++){
        data[i] = values[init + i];
    }
}

RealVector::RealVector(int init, int size, const RealVector &orig) : size_(size){
    data = new double[size];
    
    for (int i = 0; i < size; i++){
        data[i] = orig(init + i);
    }
}

RealVector::RealVector(const RealVector &orig) : size_(orig.size_) {
    data = new double[orig.size_];
    for (int i = 0; i < orig.size_; i++) data[i] = orig.data[i];
}

RealVector::RealVector(const RealVector *orig) : size_(orig->size_) {
    data = new double[orig->size_];
    for (int i = 0; i < orig->size_; i++) data[i] = orig->data[i];
}

RealVector::~RealVector(void) {
    if (data != NULL) delete [] data;
}

int RealVector::size(void)const {
    return size_;
}

void RealVector::resize(int n) {
    if (size_ == 0) {
        data = new double[n];
        for (int i = 0; i < n; i++) data[i] = 0.0;
    } else {
        double *temp0 = new double[n];

        for (int i = 0; i < min(size_, n); i++) temp0[i] = data[i];

        delete [] data;
        data = temp0;
    }

    size_ = n;

    return;
}

// Access to individual elements
//
double & RealVector::component(int n) {
    if (n < size_) return data[n];
}

const double & RealVector::component(int n)const {
    if (n < size_) return data[n];
}

// Access to data pointer
//
double * RealVector::components(void) {
    return data;
}

const double * RealVector::components(void) const {
    return data;
}

// Assignment

RealVector RealVector::operator=(const RealVector &orig) {
    // Avoid self-assignment
    if (this != &orig) {
        int n = orig.size_;
        resize(n);
        for (int i = 0; i < n; i++) data[i] = orig.data[i];
    }

    return *this;
}

bool RealVector::operator==(const RealVector &other){
    if (size_ != other.size_) return false;

    for (int i = 0; i < other.size_; i++) {
        if (data[i] != other.data[i]) return false;
    }

    return true;
}

// Cast operator
RealVector::operator double *(void) {
    return data;
}

double RealVector::operator()(int comp) const {
    return data[comp];
}


double & RealVector::operator()(int comp) {
    return data[comp];
}

double RealVector::operator[](int comp) const {
    return data[comp];
}


double & RealVector::operator[](int comp) {
    return data[comp];
}

RealVector RealVector::zeroes(int m){
    RealVector z(m);
    for (int i = 0; i < m; i++) z(i) = 0.0;

    return z;
}

// Output to stream

std::ostream & operator<<(std::ostream &out, const RealVector &r) {
    out << "(";
    for (int i = 0; i < r.size_; i++) {
        out << r.data[i];
        if (i != r.size_ - 1) out << ", ";
    }
    out << ")";

    return out;
}

// Multiplication by a scalar

RealVector operator*(const RealVector &r, double alpha) {
    RealVector temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] *= alpha;

    return temp;
}

RealVector operator*(double alpha, const RealVector &r) {
    return r*alpha;
}

// Division by a scalar
RealVector operator/(const RealVector &r, double alpha) {
    return r*(1.0/alpha);
}

// Sum with a scalar

RealVector operator+(const RealVector &r, double alpha) {
    RealVector temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] += alpha;

    return temp;
}

RealVector operator+(double alpha, const RealVector &r) {
    return r + alpha;
}

// Subtraction of/from a scalar

RealVector operator-(const RealVector &r, double alpha) {
    RealVector temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] -= alpha;

    return temp;
}

RealVector operator-(double alpha, const RealVector &r) {
    return r - alpha;
}

// Negation

RealVector operator-(const RealVector &r) {
    RealVector temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] *= -1.0;

    return temp;
}

// Sum of two RealVectors

RealVector operator+(const RealVector &x, const RealVector &y) {
    RealVector temp(x);
    for (int i = 0; i < x.size_; i++) temp.data[i] += y.data[i];

    return temp;
}

// Subtraction of two RealVectors

RealVector operator-(const RealVector &x, const RealVector &y) {
    RealVector temp(x);
    for (int i = 0; i < x.size_; i++) temp.data[i] -= y.data[i];

    return temp;
}

// Euclidean norm of a RealVector
double norm(const RealVector &x){
    return sqrt(x*x);
}

double norm2_squared(const RealVector &x){
    return x*x;
}

// Norm in L1.
double norm_L1(const RealVector &x){
    double d = 0.0;

    for (int i = 0; i < x.size(); i++) d += std::abs(x.data[i]);

    return d;
}

// Norm in L \infty.
double norm_inf(const RealVector &x){
    double d = 0.0;

    for (int i = 0; i < x.size(); i++) d = std::max(d, std::abs(x.data[i]));

    return d;
}

// Normalize a RealVector
RealVector normalize(RealVector &x){
    double inv_norm_x = 1.0/norm(x);

    x = x*inv_norm_x;
    
    return x;
}

// Inner product of two RealVectors
double operator*(const RealVector &x, const RealVector &y){
    double p = 0.0;

    for (int i = 0; i < x.size(); i++) p += x(i)*y(i);

    return p;
}

// Vector product of two 3D RealVectors
RealVector vector_product(const RealVector &x, const RealVector &y){
    RealVector v(3);

    v(0) = x(1)*y(2) - x(2)*y(1);
    v(1) = x(2)*y(0) - x(0)*y(2);
    v(2) = x(0)*y(1) - x(1)*y(0);

    return v;
}

// Solve the system of linear equations A*x = b
int solve(const DoubleMatrix &A, const RealVector &b, RealVector &x){
    int n = A.rows();

    DoubleMatrix bb(n, 1), xx(n, 1);
    for (int i = 0; i < n; i++) bb(i, 0) = b(i);

    int info = solve(A, bb, xx);

    x.resize(n);
    for (int i = 0; i < n; i++) x(i) = xx(i, 0);

    if (info == 0) return REALVECTOR_SOLVE_LINEAR_SYSTEM_OK;
    else           return REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR;
}

// Multiplication of a DoubleMatrix by a column RealVector
RealVector operator*(const DoubleMatrix &A, const RealVector &x){
    int m = A.rows(), n = A.cols();

    RealVector b(m);

    for (int i = 0; i < m; i++){
        b(i) = 0.0;
        for (int j = 0; j < n; j++) b(i) += A(i, j)*x(j);
    }

    return b;
}

// Multiplication of a row RealVector by a DoubleMatrix
RealVector operator*(const RealVector &x, const DoubleMatrix &A){
    int m = A.rows(), n = A.cols();

    RealVector b(n);

    for (int i = 0; i < n; i++){
        b(i) = 0.0;
        for (int j = 0; j < m; j++) b(i) += x(j)*A(j, i);
    }

    return b;
}

RealVector matrix_row(const DoubleMatrix &m, int r){
    int c = m.cols();

    RealVector v(c);
    for (int i = 0; i < c; i++) v(i) = m(r, i);

    return v;
}

RealVector matrix_column(const DoubleMatrix &m, int c){
    int r = m.rows();

    RealVector v(r);
    for (int i = 0; i < r; i++) v(i) = m(i, c);

    return v;
}

double vector_product_2D(const RealVector &c, const RealVector &p, const RealVector &q){
    return (p(0) - c(0))*(q(1) - c(1)) - (p(1) - c(1))*(q(0) - c(0));
}

// http://bbs.dartmouth.edu/~fangq/MATH/download/source/Determining%20if%20a%20point%20lies%20on%20the%20interior%20of%20a%20polygon.htm
//
double evaluate_line_equation(const RealVector &p0, const RealVector &p1, const RealVector &p){
    double x = p(0);
    double y = p(1);

    double x0 = p0(0);
    double y0 = p0(1);

    double x1 = p1(0);
    double y1 = p1(1);

    // val = (y - y0)*(x1 - x0) - (x - x0)*(y1 - y0)
    return (y - y0)*(x1 - x0) - (x - x0)*(y1 - y0);
}

// 
bool inside_convex_polygon(const std::vector<RealVector> &polygon, const RealVector &point){
    int n = polygon.size();

    if (n < 3) return false;

    // Obtain a first value for the equation of the line between p0 and p1 evaluated in the given point.
    //
    double val = evaluate_line_equation(polygon[0], polygon[1], point);

    bool is_inside = true;
    int pos = 1;

    // The point is inside the polygon if the line segments that form the polygon,
    // when evaluated in the given point, are all positive or negative.
    //  
    while (is_inside && pos < n){
        if (val*evaluate_line_equation(polygon[pos], polygon[(pos + 1) % n], point) < 0.0) is_inside = false;
        pos++;
    }

    return is_inside;
}

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
//
// http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
//
void convex_hull(const std::vector<RealVector> &original_polygon, std::vector<RealVector> &ch){
    std::vector<RealVector> polygon = original_polygon;

    int n = polygon.size(), k = 0;
    ch.resize(2*n);
 
    // Sort points lexicographically
    sort(polygon.begin(), polygon.end());
 
    // Build lower hull
    for (int i = 0; i < n; i++) {
        while (k >= 2 && vector_product_2D(ch[k - 2], ch[k - 1], polygon[i]) <= 0) k--;
        ch[k++] = polygon[i];
    }
 
    // Build upper hull
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && vector_product_2D(ch[k - 2], ch[k - 1], polygon[i]) <= 0) k--;
        ch[k++] = polygon[i];
    }
 
    ch.resize(k);

    return;
}

RealVector project_point_onto_line_2D(const RealVector &q, const RealVector &p0, const RealVector &p1){
    RealVector p1_minus_p0 = p1 - p0;

    DoubleMatrix P(2, 2);
    P(0, 0) = p1_minus_p0(0);
    P(0, 1) = p1_minus_p0(1);
    P(1, 0) = p1_minus_p0(1);
    P(1, 1) = -p1_minus_p0(0);

    RealVector b(2);
    b(0) = q*p1_minus_p0;
    b(1) = -p0(1)*p1_minus_p0(0) + p0(0)*p1_minus_p0(1);

    RealVector p;
    solve(P, b, p);

    return p;
}

double distance_point_line_2D(const RealVector &q, const RealVector &p0, const RealVector &p1){
    double a, b, c;
    a = p1(1) - p0(1);
    b = p0(0) - p1(0);
    c = p1(0)*p0(1) - p0(0)*p1(1);

    return std::abs(a*q(0) + b*q(1) + c)/sqrt(a*a + b*b);
}

bool segment_segment_intersection(const RealVector &p0, const RealVector &p1, const RealVector &q0, const RealVector &q1, RealVector &r, double &alpha, double &beta){
    DoubleMatrix A(2, 2);
    for (int i = 0; i < 2; i++){
        A(i, 0) = p0(i) - p1(i);
        A(i, 1) = q1(i) - q0(i);
    }

    RealVector b = q1 - p1;

    double delta = A(0, 0)*A(1, 1) - A(0, 1)*A(1, 0);
    if (fabs(delta) < 1e-10) {
        return false;
    }

    alpha = (b(0)*A(1, 1) - b(1)*A(0, 1))/delta;
    beta  = (b(1)*A(0, 0) - b(0)*A(1, 0))/delta;

    r = .5*(alpha*p0 + (1.0 - alpha)*p1 + beta*q0 + (1.0 - beta)*q1);

    return (alpha >= 0.0 && alpha <= 1.0) && (beta >= 0.0 && beta <= 1.0);
}

bool inside_non_convex_polygon(const std::vector<RealVector> &polygon, const RealVector &point){
    int j = polygon.size() - 1;
    bool oddNodes = false;

    double x = point(0);
    double y = point(1);

    for (int i = 0; i < polygon.size(); i++) {
        if ((polygon[i](1) < y && polygon[j](1) >= y || 
             polygon[j](1) < y && polygon[i](1) >= y)
            && 
            (polygon[i](0) <= x || polygon[j](0) <= x)
           ) {
            oddNodes ^= (polygon[i](0) + (y - polygon[i](1))/(polygon[j](1) - polygon[i](1))*(polygon[j](0) - polygon[i](0)) < x); 
        }
        j = i;
    }

    return oddNodes;
}

double sum(const RealVector &v){
    double s = 0.0;

    for (int i = 0; i < v.size(); i++) s += v(i);

    return s;
}


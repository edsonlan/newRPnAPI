#include "PointND.h"
#include "BoxND.h"

inline int PointND::min(int x, int y){
    return (x < y) ? x : y; 
}

PointND::PointND(void){
    data  = NULL;
    size_ = 0;
}

PointND::PointND(int n) : size_(n){
    data = new double[n];
    for (int i = 0; i < n; i++) data[i] = 0.0;
}

PointND::PointND(const PointND &orig) : size_(orig.size_){
    data = new double[orig.size_];
    for (int i = 0; i < orig.size_; i++) data[i] = orig.data[i];
}

PointND::PointND(const PointND *orig) : size_(orig->size_){
    data = new double[orig->size_];
    for (int i = 0; i < orig->size_; i++) data[i] = orig->data[i];
}

PointND::~PointND(void){
    if (data != NULL) delete [] data;
}

int PointND::size(void) const {
    return size_;
}

void PointND::resize(int n){
    if (size_ == 0){
        data = new double[n];
        for (int i = 0; i < n; i++) data[i] = 0.0;
    }
    else {
        double *temp0 = new double[n];

        for (int i = 0; i < min(size_, n); i++) temp0[i] = data[i];

        delete [] data;
        data = temp0;
    }

    size_ = n;

    return;
}

// Access to individual elements
double & PointND::component(int n){
    if (n >= 0 && n < size_) return data[n];
}

double PointND::component(int n) const {
    if (n >= 0 && n < size_) return data[n];
}

double & PointND::operator()(int n){
    return component(n);
}

double PointND::operator()(int n) const {
    return component(n);
}

// Access to data pointer
double * PointND::components(void){
    return data;
}

double * PointND::components(void) const {
    return data;
}

// Assignment
PointND PointND::operator=(const PointND &orig){
    // Avoid self-assignment
    if (this != &orig){
        int n = orig.size_;
        resize(n);
        for (int i = 0; i < n; i++) data[i] = orig.data[i];
    }

    return *this;
}

// Output to stream
std::ostream & operator<<(std::ostream &out, const PointND &r){
    out << "(";
    for (int i = 0; i < r.size_; i++){
        out << r.data[i];
        if (i != r.size_ - 1) out << ", ";
    }
    out << ")";

    return out;
}

// Multiplication by a scalar
PointND operator*(const PointND &r, double alpha){
    PointND temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] *= alpha;

    return temp;
}

PointND operator*(double alpha, const PointND &r){
    return r*alpha;
}

// Division by a scalar
PointND operator/(const PointND &r, double alpha){
    PointND temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] /= alpha;

    return temp;
}

// Sum with a scalar
PointND operator+(const PointND &r, double alpha){
    PointND temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] += alpha;

    return temp;
}

PointND operator+(double alpha, const PointND &r){
    return r + alpha;
}

// Subtraction of/from a scalar
PointND operator-(const PointND &r, double alpha){
    PointND temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] -= alpha;

    return temp;
}


PointND operator-(double alpha, const PointND &r){
    return r - alpha;
}

// Negation
PointND operator-(const PointND &r){
    PointND temp(r);
    for (int i = 0; i < r.size_; i++) temp.data[i] *= -1.0;

    return temp;
}

// Sum of two PointNDs
PointND operator+(const PointND &x, const PointND &y){
    PointND temp(x);
    for (int i = 0; i < x.size_; i++) temp.data[i] += y.data[i];

    return temp;
}

// Subtraction of two PointNDs
PointND operator-(const PointND &x, const PointND &y){
    PointND temp(x);
    for (int i = 0; i < x.size_; i++) temp.data[i] -= y.data[i];

    return temp;
}

bool PointND::intersect(const BoxND &b) const {
    bool intersection = true;
    int i = 0;

    while (intersection && i < size_){
        intersection = (data[i] >= b.pmin.data[i]) && (data[i] <= b.pmax.data[i]);
        i++;
    }

    return intersection;
}


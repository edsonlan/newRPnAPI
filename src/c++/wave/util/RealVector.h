#ifndef _REALVECTOR_
#define _REALVECTOR_

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "DoubleMatrix.h"

#ifndef REALVECTOR_SOLVE_LINEAR_SYSTEM_OK
#define REALVECTOR_SOLVE_LINEAR_SYSTEM_OK 0
#endif

#ifndef REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR
#define REALVECTOR_SOLVE_LINEAR_SYSTEM_ERROR (-1)
#endif

class RealVector {
private:
protected:
    double *data;
    int size_;

    int min(int x, int y);
public:
    RealVector(void);
    RealVector(int n);
    RealVector(int n, const double *);
    RealVector(int init, int size, double *);
    RealVector(int init, int size, const RealVector &orig);
    RealVector(const RealVector &orig);
    RealVector(const RealVector *orig);

    ~RealVector(void);

    int size(void)const;
    void resize(int n);

    // Access to individual elements
    double & component(int n);
    const double & component(int n) const;

    double operator()(int comp) const;
    double & operator()(int comp);

    double operator[](int comp) const;
    double & operator[](int comp);

    // Cast operator
    operator double *(void);

    // Access to data pointer
    double * components(void);
    const double * components(void) const;

    // Assignment
    RealVector operator=(const RealVector &orig);
    bool operator==(const RealVector &other);

    // Return a vector of zeroes
    static RealVector zeroes(int m);

    // Output to stream
    friend std::ostream & operator<<(std::ostream &out, const RealVector &r);

    // Multiplication by a scalar
    friend RealVector operator*(const RealVector &r, double alpha);
    friend RealVector operator*(double alpha, const RealVector &r);

    // Division by a scalar
    friend RealVector operator/(const RealVector &r, double alpha);

    // Sum with a scalar
    friend RealVector operator+(const RealVector &r, double alpha);
    friend RealVector operator+(double alpha, const RealVector &r);

    // Subtraction of/from a scalar
    friend RealVector operator-(const RealVector &r, double alpha);
    friend RealVector operator-(double alpha, const RealVector &r);

    // Negation
    friend RealVector operator-(const RealVector &r);

    // Sum of two RealVectors
    friend RealVector operator+(const RealVector &x, const RealVector &y);

    // Subtraction of two RealVectors
    friend RealVector operator-(const RealVector &x, const RealVector &y);
 
    // Euclidean norm of a RealVector
    friend double norm(const RealVector &x);
    friend double norm2_squared(const RealVector &x);

    // Norm in L1.
    friend double norm_L1(const RealVector &x);

    // Norm in L \infty.
    friend double norm_inf(const RealVector &x);

    // Normalize a RealVector
    friend RealVector normalize(RealVector &x);

    // Inner product of two RealVectors
    friend double operator*(const RealVector &x, const RealVector &y);

    // Vector product of two 3D RealVectors
    friend RealVector vector_product(const RealVector &x, const RealVector &y);

    // Solve the system of linear equations A*x = b
    friend int solve(const DoubleMatrix &A, const RealVector &b, RealVector &x);

    // Multiplication of a DoubleMatrix by a column RealVector
    friend RealVector operator*(const DoubleMatrix &A, const RealVector &x);

    // Multiplication of a row RealVector by a DoubleMatrix
    friend RealVector operator*(const RealVector &x, const DoubleMatrix &A);

    RealVector& operator+=(const RealVector &v){
        for (int i = 0; i < size(); i++) component(i) += v(i);

        return *this;
    }

    RealVector& operator-=(const RealVector &v){
        for (int i = 0; i < size(); i++) component(i) -= v(i);

        return *this;
    }

    RealVector& operator+=(double v){
        for (int i = 0; i < size(); i++) component(i) += v;

        return *this;
    }

    RealVector& operator-=(double v){
        for (int i = 0; i < size(); i++) component(i) -= v;

        return *this;
    }

    RealVector& operator*=(double v){
        for (int i = 0; i < size(); i++) component(i) *= v;

        return *this;
    }

    RealVector& operator/=(double v){
        for (int i = 0; i < size(); i++) component(i) /= v;

        return *this;
    }

    // The following methods are used when computing a convex hull in 2D and if
    // a point is inside a convex hull. The extension will use them.

    // To sort lexicographically
    bool operator<(const RealVector &p) const {
        return std::lexicographical_compare(data, data + size_, p.data, p.data + p.size_);
    }

    // Vector product in 2D of cp and cq.
    //
    friend double vector_product_2D(const RealVector &c, const RealVector &p, const RealVector &q);

    // http://bbs.dartmouth.edu/~fangq/MATH/download/source/Determining%20if%20a%20point%20lies%20on%20the%20interior%20of%20a%20polygon.htm
    //
    friend double evaluate_line_equation(const RealVector &p0, const RealVector &p1, const RealVector &p);

    // 
    friend bool inside_convex_polygon(const std::vector<RealVector> &polygon, const RealVector &point);

    // Returns a list of points on the convex hull in counter-clockwise order.
    // Note: the last point in the returned list is the same as the first one.
    //
    // http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
    //
    friend void convex_hull(const std::vector<RealVector> &original_polygon, std::vector<RealVector> &ch);

    // Based on
    //
    //     http://cs.nyu.edu/~yap/classes/visual/03s/hw/h2/math.pdf
    //
    friend RealVector project_point_onto_line_2D(const RealVector &q, const RealVector &p0, const RealVector &p1);

    // Based on
    //
    //     http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    //
    friend double distance_point_line_2D(const RealVector &q, const RealVector &p0, const RealVector &p1);

    // r = alpha*p0 + (1.0 - alpha)*p1 = beta*q0 + (1.0 - beta)*q1.
    //
    friend bool segment_segment_intersection(const RealVector &p0, const RealVector &p1, const RealVector &q0, const RealVector &q1, RealVector &r, double &alpha, double &beta);

    // Point inside a non-convex polygon. Based on:
    //
    //    http://alienryderflex.com/polygon/
    //
    friend bool inside_non_convex_polygon(const std::vector<RealVector> &polygon, const RealVector &point); 

    // Sum of all the elements of a vector.
    //
    friend double sum(const RealVector &v);
};

// Extract rows and columns of a DoubleMatrix and return them as RealVectors.
RealVector matrix_row(const DoubleMatrix &m, int r);
RealVector matrix_column(const DoubleMatrix &m, int c);

#endif // _REALVECTOR_


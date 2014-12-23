#ifndef _POINTND_
#define _POINTND_

#include <iostream>
#include <fstream>

class BoxND;

class PointND {
    private:
    protected:
        double *data;
        int size_;

        int min(int x, int y);
    public:
        PointND(void);
        PointND(int n);
        PointND(const PointND &orig);
        PointND(const PointND *orig);

        ~PointND(void);

        int size(void) const;
        void resize(int n);

        // Access to individual elements
        double       & component(int n);
        double         component(int n) const;
        double       & operator()(int n);
        double         operator()(int n) const;

        // Access to data pointer
        double * components(void);
        double * components(void) const;

        // Assignment
        PointND operator=(const PointND &orig);

        // Output to stream
        friend std::ostream & operator<<(std::ostream &out, const PointND &r);

        // Multiplication by a scalar
        friend PointND operator*(const PointND &r, double alpha);
        friend PointND operator*(double alpha, const PointND &r);

        // Division by a scalar
        friend PointND operator/(const PointND &r, double alpha);

        // Sum with a scalar
        friend PointND operator+(const PointND &r, double alpha);
        friend PointND operator+(double alpha, const PointND &r);

        // Subtraction of/from a scalar
        friend PointND operator-(const PointND &r, double alpha);
        friend PointND operator-(double alpha, const PointND &r);

        // Negation
        friend PointND operator-(const PointND &r);

        // Sum of two PointNDs
        friend PointND operator+(const PointND &x, const PointND &y);

        // Subtraction of two PointNDs
        friend PointND operator-(const PointND &x, const PointND &y);

        // Intersection with a box
        bool intersect(const BoxND &b) const ;
};

#endif // _POINTND_


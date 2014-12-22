#ifndef _UTILITIES_
#define _UTILITIES_

#include <limits>
#include <vector>
#include "RealVector.h"
#include "WaveCurve.h"

#define BISECTIONONSEGMENT_OK    0
#define BISECTIONONSEGMENT_ERROR 1

#define UTILITIES_BISECTION_OK    0
#define UTILITIES_BISECTION_ERROR 1

#define BHASKARA_COMPLEX_ROOTS       (-1)
#define BHASKARA_DOUBLE_ROOTS          0
#define BHASKARA_TWO_DIFFERENT_ROOTS   1

class Utilities {
    private:
    protected:
        static bool validate_alpha(double alpha){
            return (alpha >= 0.0 && alpha <= 1.0);
        }
    public:
        // Convert a vector of RealVectors into a vector of segments.
        //
        static void points_to_segments(std::vector<RealVector> &v);

        // Bisection to find the zero of function_for_bisection, given a generic object obj.
        //
        static int bisection_on_segment(void* obj, double (*function_for_bisection)(void*, const RealVector&), const RealVector &p, const RealVector &q, RealVector &r);

        static void alpha_convex_combination_projection(const RealVector &p0, const RealVector &p1, const RealVector &p, double &alpha, RealVector &proj);

        static void pick_point_from_continuous_curve(const std::vector<RealVector> &curve, const RealVector &p, 
                                                     RealVector &closest_point);
        static void pick_point_from_continuous_curve(const std::vector<RealVector> &curve, const RealVector &p, 
                                                     RealVector &closest_point, int &closest_segment_index);

        static void pick_point_from_segmented_curve(const std::vector<RealVector> &curve, const RealVector &p, RealVector &closest_point);
        static void pick_point_from_segmented_curve(const std::vector<RealVector> &curve, const RealVector &p, RealVector &closest_point, int &index_p0, int &index_p1);

        static void pick_point_from_wavecurve(const WaveCurve &wavecurve, const RealVector &p, 
                                              int &curve_index, int &segment_index_in_curve, RealVector &closest_point, double &speed);

        // Created specifically for the Injection_to_Side case (Freddie's).
        //
        static void pick_last_point_from_wavecurve(const WaveCurve &wavecurve, const RealVector &p, 
                                                   int &curve_index, int &segment_index_in_curve, RealVector &closest_point, double &speed);

        static int find_point_on_level_curve(void *obj, double (*function_for_bisection)(void*, const RealVector &), const RealVector &p, const RealVector &q0, RealVector &point_on_level_curve);

        // Subdivide a segment with vertices p and q in a Curve with (n + 1) vertices.
        // 
        static void regularly_sampled_segment(const RealVector &p, const RealVector &q, int n, Curve &curve);

        // Find the roots of the polynomial x*x + b*x + c.
        //
        static int Bhaskara(double b, double c, double &x1, double &x2);
        static int Bhaskara(double a, double b, double c, double &x1, double &x2);

        template <typename T> static int sgn(T val) {
            return (T(0) < val) - (val < T(0));
        }
};

#endif // _UTILITIES_


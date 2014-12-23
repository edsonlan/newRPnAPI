#ifndef _CONTOUR2X2_METHOD_
#define _CONTOUR2X2_METHOD_

#include "HyperCube.h"
#include "RealSegment.h"
#include "ThreeImplicitFunctions.h"

#include <iostream>
#include <vector>

#include "Boundary.h"
#include "GridValues.h"

//#ifndef INVALID_FUNCTION_ON_VERTICES
//#define INVALID_FUNCTION_ON_VERTICES 0
//#endif

//#ifndef VALID_FUNCTION_ON_VERTICES
//#define VALID_FUNCTION_ON_VERTICES 1
//#endif

//#ifndef NO_ZERO_ON_CUBE
//#define NO_ZERO_ON_CUBE 2
//#endif

class Contour2x2_Method {
    private:
        int dimension;

        static bool is_first;

        static HyperCube hc;
        static int hn;
        static int hm;
        static int nsface_, nface_, nsoln_, nedges_;
        static int dims_;
        static int dime_;
        static int dimf_;
        static int ncvert_;
        static int dncv;  // TODO, ver si no es ncvert.
        static int nsimp_;

        static int numberOfCombinations;

        static int *storn_;
        static int *storm_;

        static Matrix<double> cvert_;
        static Matrix<double> vert;
        static Matrix<int>    bsvert_;
        static Matrix<int>    perm_;
        static Matrix<int>    comb_; 
        static Matrix<double> foncub;
        static Matrix<int>    fnbr_;

        static Matrix<int>    facptr_;
        static Matrix<int>    face_;

        static Matrix<double> cpp_sol;
        static Matrix<int>    solptr_;
        static Matrix<int>    cpp_edges_;
        static Matrix<int>    smpedg_;

        static int *exstfc;
        static int *sptr_;

        static int tsimp;
        static int tface;

        static int *index;
        static int *usevrt; // TODO, ver si es necesario, sino, matar [dncv]

        static void allocate_arrays(void);
    protected:
        static double ul0, vl0, ur0, vr0;
        static double dul, dvl, dur, dvr;
        static double dumax, dvmax;

        static bool filhcub4(ThreeImplicitFunctions *timpf,
                             int ir, int jr, int *index, double *foncub);

        static void filedg4(Matrix<double> &sol_, Matrix<int> &edges_, int nedges_, 
                            int il, int jl, int ir, int jr,
                            std::vector<RealVector> &left_vrs, std::vector<RealVector> &right_vrs);

		static bool left_right_adjacency(int il, int jl, int ir, int jr);
    public:
        static void deallocate_arrays(void);

        static void curve2x2(ThreeImplicitFunctions *timpf,
                             std::vector<RealVector> &left_vrs, 
                             std::vector<RealVector> &right_vrs);
};

#endif // _CONTOUR2X2_METHOD_


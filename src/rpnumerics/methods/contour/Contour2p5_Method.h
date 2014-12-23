#ifndef _CONTOUR2P5_METHOD_
#define _CONTOUR2P5_METHOD_

#include "HyperCube.h"
//#include "RealSegment.h"
#include "TwoImplicitFunctions.h"

#include <iostream>
#include <vector>

#include "Boundary.h"
#include "GridValues.h"

#ifndef INVALID_FUNCTION_ON_VERTICES
#define INVALID_FUNCTION_ON_VERTICES 0
#endif

#ifndef VALID_FUNCTION_ON_VERTICES
#define VALID_FUNCTION_ON_VERTICES 1
#endif

#ifndef NO_ZERO_ON_CUBE
#define NO_ZERO_ON_CUBE 2
#endif

class Contour2p5_Method {
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
        static void curve2p5(TwoImplicitFunctions *timpf, int current_segment_index, 
                             std::vector<RealVector> &curve_vrs, 
                             std::vector<RealVector> &domain_vrs);
                             
        static void filedg3(Matrix<double> &sol_, Matrix<int> &edges_,
                            int nedges_, const RealVector &grid_resolution,
                            const std::vector<RealVector> &segment,
                            const RealVector &uv,
                            std::vector<RealVector> &curve_segments, std::vector<RealVector> &domain_segments);
                                
        static int function_on_cube(TwoImplicitFunctions *timpf, Matrix<double> &foncub, int ir, int jr);  
              
    public:
        static void deallocate_arrays(void);

        static void contour2p5(TwoImplicitFunctions *timpf, std::vector<RealVector> &curve_vrs, std::vector<RealVector> &domain_vrs);
        
        static void contour2p5_for_curve(TwoImplicitFunctions *timpf, std::vector<RealVector> &curve_vrs, std::vector<RealVector> &domain_vrs);
};

#endif // _CONTOUR2P5_METHOD_


/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) HugoniotContourMethod.h
 */

#ifndef _HugoniotContourMethod_H
#define _HugoniotContourMethod_H

#define CONTINUATION_METHOD 1
#define SEGMENTATION_METHOD 0

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "HyperCube.h"
#include "RealSegment.h"
#include "HugoniotFunctionClass.h"

#include <iostream>
#include <vector>
#include <deque>

#include "ImplicitFunction.h"
#include "Newton_Improvement.h"
#include "Boundary.h"
#include "GridValues.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

class ContourMethod {
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
    static int nsimp_;

    static int numberOfCombinations;

    static int *storn_;
    static int *storm_;

    static double *cvert_;
    static double *vert;
    static int *bsvert_;
    static int *perm_;
    static int *comb_;
    static int *fnbr_;

    static int *facptr_;
    static int *face_;

    static double *sol_;
    static int *solptr_;
    static int *edges_;
    static int *smpedg_;

    static int *exstfc;
    static int *sptr_;
    static int * gamb;

    static int tsimp;
    static int tface;

    static void allocate_arrays(void);

    // Auxiliary grid for continuation
    //
    static Matrix<int> number_chains;
    static Matrix<std::vector <std::vector <int> > > chain_edges;
    static Matrix<std::vector <std::vector <RealVector> > > chains;
//    static Matrix<bool> iplus, iminus, jplus, jminus;

    static std::vector < std::vector <int> > chain_list;

    static int topological_sort(int i, int j);


public:
    static void deallocate_arrays(void);

protected:
public:
    static int contour2d(ImplicitFunction *impf, std::vector<RealVector> &vrs);

    static int contour2d(ImplicitFunction *impf,
            std::vector< std::deque <RealVector> > &curves,
            std::vector <bool> &is_circular);

    static int contour2d(ImplicitFunction *impf, std::vector<RealVector> &vrs,
            std::vector< std::deque <RealVector> > &curves,
            std::vector <bool> &is_circular,
            const int method);

    static int contour2d(ImplicitFunction *impf, int min_row, int max_row, int min_col, int max_col,
                         std::vector<RealVector> &vrs,
                         std::vector< std::deque <RealVector> > &curves,
                         std::vector <bool> &is_circular,
                         const int method);
};

#endif //! _HugoniotContourMethod_H

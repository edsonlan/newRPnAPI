#ifndef _HugoniotContourMethodPure_H
#define _HugoniotContourMethodPure_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "HyperCube.h"
#include "RealSegment.h"
#include "HugoniotFunctionClass.h"

#include <iostream>
#include <vector>

#include "ImplicitFunction.h"
#include "Newton_Improvement.h"
#include "Boundary.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

// Forward declaration:

//class Inflection_Curve;

class ContourMethodPure {
    private:
        int dimension;

//        static int inpdom(double *u);

        HugoniotFunctionClass *hugoniot;	 // TODO: Edson, voce deve comentar essa linha (e no .cc)

        //
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

        static int tsimp;
        static int tface;

//<<<<<<< HEAD
    //void classifiedCurve(const RealVector &, int, int, vector<HugoniotPolyLine> &);
    //void unclassifiedCurve(const RealVector &, vector<HugoniotPolyLine> &);
//=======
        static void allocate_arrays(void);


    public:
        static void deallocate_arrays(void);

    protected:
    public:
        ContourMethodPure(HugoniotFunctionClass *h); // TODO: Edson, voce deve comentar essa linha (e no .cc)
        ~ContourMethodPure();			 // TODO: Edson, voce deve comentar essa linha (e no .cc)
    
        //int curv2d(int sn, int seglim, double f, double *rect, int *res, int ifirst) ;
        //int curv2d(int sn, int seglim, double f, double *rect, int *res, int ifirst, std::vector<RealSegment> &vrs) ;
        int curv2d(int sn, int seglim, double f, double *rect, int *res, int ifirst, std::vector<RealVector> &vrs);

//        static int vect2d(Inflection_Curve *functioninsquare, int &sn, int &seglim, double f, double *rect, int *res, int ifirst, std::vector<RealVector> &vrs);
        static int contour2d(ImplicitFunction *impf, Boundary *boundary, double *rect, int *res, std::vector<RealVector> &vrs);
};

#endif //! _HugoniotContourMethodPure_H


#ifndef _CURVE_
#define _CURVE_

#include <vector>
#include "RealVector.h"
#include "ReferencePoint.h"
#include "HyperOctree.h"
#include "GridValues.h"

#define RAREFACTION_CURVE 1
#define COMPOSITE_CURVE   2
#define SHOCK_CURVE       3

#define SPEED_INCREASE   10
#define SPEED_DECREASE   11

#define ALL_SIDES              20
#define ALL_SIDES_EXCEPT_LEFT  21
#define ALL_SIDES_EXCEPT_RIGHT 22
#define ALL_SIDES_EXCEPT_UP    23
#define ALL_SIDES_EXCEPT_DOWN  24

class Curve {
    private:
    protected:
        bool segment_intersects_box(const RealVector &seg0, const RealVector &seg1, const RealVector &b0, const RealVector &b1, int side_to_check, int &side);
        bool segment_intersects_box(const RealVector &seg0, const RealVector &seg1, const RealVector &box0, const RealVector &box1);
    public:
        // RAREFACTION_CURVE or COMPOSITE_CURVE or SHOCK_CURVE.
        //
        int type;

        // From whence this curve came.
        //
        int back_curve_index;

        // These lines below should be moved to WaveCurve (this impacts RarefactionCurve, CompositeCurve and maybe ShockCurve).
        int family;

        int increase;

        ReferencePoint reference_point;

        //
        //     this->curve[i] is related to wavecurve[back_curve_index][back_pointer[i]] 
        //     
        // More precisely:
        // 
        //     this->curve[i] is related to wavecurve[back_curve_pointer[i]][back_pointer[i]] 
        //
        std::vector<RealVector> curve;
        std::vector<double>     xi;
        std::vector<int>        back_pointer;
        std::vector<int>        back_curve_pointer;

        std::vector<int> explicit_bifurcation_transition_index;

        std::vector<int>        extension_of_whom;

        std::vector<bool> compatible;

        RealVector last_point;
        RealVector final_direction;

        // Why this curve stopped.
        //
        int reason_to_stop;

        // The following members are not to be filled when computing an integral curve:
        //
        std::vector<double>     speed;
        
        std::vector<RealVector> eigenvalues;

        void init(const Curve &orig);

        Curve();
        Curve(const Curve &orig);
        Curve(const Curve *orig);

        virtual ~Curve();

        Curve operator=(const Curve &orig);

        void clear();

        void disable_adjacent_cells(GridValues *g, int row, int col);
        void disable_intersecting_cells(GridValues *g);
        void disable_intersecting_cells(const RealVector &p, const RealVector &q, GridValues *g);
};

#endif // _CURVE_


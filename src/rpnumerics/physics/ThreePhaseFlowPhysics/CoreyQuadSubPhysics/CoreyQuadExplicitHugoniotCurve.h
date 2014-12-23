#ifndef _COREYQUADEXPLICITHUGONIOTCURVE_
#define _COREYQUADEXPLICITHUGONIOTCURVE_

//#include "CoreyQuadHugoniotCurve.h"
#include "HugoniotCurve.h"
#include "ParametricPlot.h"
#include "Utilities.h"

#define COREYQUAD_GENERIC_POINT 0
#define COREYQUAD_G_VERTEX      1
#define COREYQUAD_W_VERTEX      2
#define COREYQUAD_O_VERTEX      3
#define COREYQUAD_GW_SIDE       4
#define COREYQUAD_WO_SIDE       5
#define COREYQUAD_GO_SIDE       6
#define COREYQUAD_G_BIFURCATION 7
#define COREYQUAD_W_BIFURCATION 8
#define COREYQUAD_O_BIFURCATION 9
#define COREYQUAD_UMBILIC_POINT 10
#define COREYQUAD_E_POINT       11
#define COREYQUAD_B_POINT       12
#define COREYQUAD_D_POINT       13

class CoreyQuadSubPhysics;

class CoreyQuadExplicitHugoniotCurve : public HugoniotCurve {
    private:
    protected:
        static double sign(double v){
            return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0.0);
        }

        std::vector<int> side_opposite_vertex;

        double muw, mug, muo;

        // This pointer will take the address of one of the methods below.
        //
        RealVector (*function_to_use)(void *cqeh, double phi);

        static RealVector generic(void *cqeh, double phi);

//        static RealVector G_vertex(void *cqeh, double phi);
//        static RealVector O_vertex(void *cqeh, double phi);
//        static RealVector W_vertex(void *cqeh, double phi);

        static RealVector G_bif(void *cqeh, double phi);
        static RealVector O_bif(void *cqeh, double phi);
        static RealVector W_bif(void *cqeh, double phi);

        static RealVector GO_side(void *cqeh, double phi);
        static RealVector GW_side(void *cqeh, double phi);
        static RealVector WO_side(void *cqeh, double phi);

//        static RealVector umbilic(void *cqeh, double phi);

        RealVector compute_umbilic_point(const RealVector &mu);

        void vertex_and_side(int side_opposite_vertex, const RealVector &mu, RealVector &vertex, RealVector &point_on_side);

        CoreyQuadSubPhysics *corey;
    public:
        CoreyQuadExplicitHugoniotCurve(CoreyQuadSubPhysics *c);
        virtual ~CoreyQuadExplicitHugoniotCurve();

        void curve(const ReferencePoint &ref, int type, std::vector<Curve> &c);

//        void curve(const ReferencePoint &ref, int type, std::vector<HugoniotPolyLine> &classified_curve){
//            CoreyQuadHugoniotCurve::curve(ref, type, classified_curve);

//            return;
//        }

        // To be used by ParametricPlot.
        //
//        static RealVector f(void *obj, double phi);
        static bool f_asymptote(void *obj, const RealVector &p, const RealVector &q);

        void list_of_reference_points(std::vector<int> &type, std::vector<std::string> &name) const {
            type.clear();

            type.push_back(COREYQUAD_GENERIC_POINT);
            type.push_back(COREYQUAD_G_VERTEX);
            type.push_back(COREYQUAD_W_VERTEX);
            type.push_back(COREYQUAD_O_VERTEX);
            type.push_back(COREYQUAD_GW_SIDE);
            type.push_back(COREYQUAD_WO_SIDE);
            type.push_back(COREYQUAD_GO_SIDE);
            type.push_back(COREYQUAD_G_BIFURCATION);
            type.push_back(COREYQUAD_W_BIFURCATION);
            type.push_back(COREYQUAD_O_BIFURCATION);
            type.push_back(COREYQUAD_UMBILIC_POINT);
            type.push_back(COREYQUAD_E_POINT);
            type.push_back(COREYQUAD_B_POINT);
            type.push_back(COREYQUAD_D_POINT);

            name.clear();

            name.push_back(std::string("Generic point"));
            name.push_back(std::string("G vertex"));
            name.push_back(std::string("W vertex"));
            name.push_back(std::string("O vertex"));
            name.push_back(std::string("GW side"));
            name.push_back(std::string("WO side"));
            name.push_back(std::string("GO side"));
            name.push_back(std::string("G bifurcation"));
            name.push_back(std::string("W bifurcation"));
            name.push_back(std::string("O bifurcation"));
            name.push_back(std::string("Umbilic point"));
            name.push_back(std::string("E")); 
            name.push_back(std::string("B"));   
            name.push_back(std::string("D"));

            return;
        }
};

#endif // _COREYQUADEXPLICITHUGONIOTCURVE_


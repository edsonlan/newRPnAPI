#ifndef _COMPOSITECURVE_
#define _COMPOSITECURVE_


#include <cmath> // For std::abs.
#include <algorithm>


#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Boundary.h"
#include "ShockCurve.h"
#include "ReferencePoint.h"

#include "Explicit_Bifurcation_Curves.h"

#include "ODE_Solver.h"
#include "Bisection.h"


// Forward declaration.
//
class WaveCurve;

#ifndef COMPOSITE_OK
#define COMPOSITE_OK    0
#endif

#ifndef COMPOSITE_ERROR
#define COMPOSITE_ERROR 1
#endif

// Reasons of error.
//
#ifndef COMPOSITE_ERROR_AT_BEGINNING_OUT_OF_BOUNDARY
#define COMPOSITE_ERROR_AT_BEGINNING_OUT_OF_BOUNDARY 2
#endif

#ifndef COMPOSITE_ERROR_AT_BEGINNING_DETERMINANT
#define COMPOSITE_ERROR_AT_BEGINNING_DETERMINANT 3
#endif

#ifndef COMPOSITE_ERROR_AT_DETERMINANT
#define COMPOSITE_ERROR_AT_DETERMINANT 4
#endif

#ifndef COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING
#define COMPOSITE_ERROR_AT_RAREFACTION_BEGINNING 5
#endif

#ifndef COMPOSITE_LAST_POINT_ERROR
#define COMPOSITE_LAST_POINT_ERROR 6
#endif

// Options.
//
#ifndef COMPOSITE_BEGINS_AT_INFLECTION
#define COMPOSITE_BEGINS_AT_INFLECTION 10
#endif

#ifndef COMPOSITE_AFTER_COMPOSITE
#define COMPOSITE_AFTER_COMPOSITE      11
#endif

// Reasons of success
#ifndef COMPOSITE_REACHED_BOUNDARY
#define COMPOSITE_REACHED_BOUNDARY       100
#endif

#ifndef COMPOSITE_REACHED_DOUBLE_CONTACT
#define COMPOSITE_REACHED_DOUBLE_CONTACT 101
#endif

#ifndef COMPOSITE_REACHED_LINE
#define COMPOSITE_REACHED_LINE           102
#endif

// Equivalent to "Rarefaction depleted."
//
#ifndef COMPOSITE_COMPLETED    
#define COMPOSITE_COMPLETED              103
#endif

#ifndef SECUNDARY_BIFURCATION_DETECTED    
#define SECUNDARY_BIFURCATION_DETECTED       200
#endif

#ifndef SECUNDARY_BIFURCATION_NOT_DETECTED    
#define SECUNDARY_BIFURCATION_NOT_DETECTED   201
#endif

#ifndef NORMALIZE_WITH_RESPECT_TO_COMPOSITE
#define NORMALIZE_WITH_RESPECT_TO_COMPOSITE 300
#endif

#ifndef NORMALIZE_WITH_RESPECT_TO_RAREFACTION
#define NORMALIZE_WITH_RESPECT_TO_RAREFACTION 301
#endif

    class alpha_index {
        public:
            double alpha;
            int index;

            alpha_index() : alpha(0.0), index(0){}

            alpha_index(double a, int i) : alpha(a), index(i){}

            ~alpha_index(){}

            friend bool operator<(const alpha_index &a, const alpha_index &b){
                return a.alpha < b.alpha;
            }

            alpha_index & operator=(const alpha_index &orig){
                if (&orig != this){
                    alpha = orig.alpha;
                    index = orig.index;
                }

                return *this;
            }
    };

class CompositeCurve {
    private:
    protected:
        const FluxFunction *flux;
        const AccumulationFunction *accum;
        const Boundary *boundary;
        ShockCurve *shock;

        int retreat;

        // For integrating as a ODE.
        int family;        
        double tolerance;
        RealVector reference_vector;
        double referencedeterminant; // If referencedeterminant
        
        // TODO: This should be replaced by a pointer to the rarefaction curve as a whole. A class RarefactionCurvePoints should be created
        //       that holds lots of information: points, lambdas, family, etc.
        double lambda_init_base_rarefaction;

        const FluxFunction         *rarflux;
        const AccumulationFunction *raraccum;
        const Boundary             *rarboundary;

        RealVector composite_field(const RealVector &final_point_pair);

        virtual void all_eigenvalues(const RealVector &p, int family, RealVector &point_eigenvalues);
        virtual void add_point_to_curve(const RealVector &p, int back, const Curve &rarcurve, Curve &curve);

        // To be used by the correction of the last point.
        //                #include "Hugoniot_Curve.h"
        double lambda_at_double_contact;
        RealVector rar_F_at_double_contact, rar_G_at_double_contact;

        // If this composite is not the first one and a previous composite 
        // reached a double contact, it may happen that this one will exhaust its
        // companion rarefaction. To compute the last point of such a composite
        // a different signal event will be used. Find the corresponding 
        // signal event method in the public interface (below).
        //
        double maxsigma;
        bool use_maxsigma;
        
        // TESTING.
        bool compute_first_determinant;
        double first_determinant;
        int (*field)(int *two_n, double *xi, double *pointpair, double *field, int *obj, double *data);
        double cmp_deltaxi;
        bool use_field_near_double_contact;

        // Explicit bifurcation transitions.
        //
        Explicit_Bifurcation_Curves *explicit_bifurcation_curve;
        int index_of_explicit_bifurcation_expression;
//        void transition_with_explicit_bifurcation(const RealVector &rarcmp_point, double initial_time, RealVector &out, double &final_time, int &index_of_explicit_bifurcation);

        int normalize_with_respect_to_whom;

       
    public:
        CompositeCurve(const AccumulationFunction *a, const FluxFunction *f, const Boundary *b, ShockCurve *s, Explicit_Bifurcation_Curves *ebc);
        virtual ~CompositeCurve();

        static int composite_field(int *two_n, double *xi, double *pointpair, double *field, int *obj, double *data);
        static int composite_field_near_double_contact(int *two_n, double *xi, double *pointpair, double *field, int *obj, double *data);

        static int composite_field_near_double_contact3D2D(int *two_n, double *xi, double *pointpair, double *field, int *obj, double* /* Not used */);

//        int curve(const FluxFunction *RarFlux, const AccumulationFunction *RarAccum, 
//                  const Boundary *RarBoundary, std::vector<RealVector> &rarcurve, std::vector<double> &lambda,
//                  const RealVector &composite_initial_point,
//                  const ODE_Solver *odesolver,
//                  double deltaxi,
//                  int where_composite_begins, int fam, 
//                  std::vector<RealVector> &newrarcurve,
//                  std::vector<RealVector> &compositecurve,
//                  RealVector &final_direction,
//                  std::vector<double> &dets,
//                  int &reason_why,
//                  int &edge);

        virtual int curve(const AccumulationFunction *RarAccum, const FluxFunction *RarFlux,
                  const Boundary *RarBoundary, 
                  const Curve &rarcurve,
//                  std::vector<RealVector> &rarcurve, std::vector<double> &lambda,
                  const RealVector &composite_initial_point,
                  int last_point_in_rarefaction,
                  const ODE_Solver *odesolver,
                  double deltaxi,
                  int where_composite_begins, int fam, 
//                  std::vector<RealVector> &newrarcurve,
//                  std::vector<RealVector> &compositecurve,
                  Curve &new_rarcurve,
                  Curve &compositecurve,
                  RealVector &final_direction,
                  int &reason_why,
                  int &edge){

            return curve(RarAccum, RarFlux, RarBoundary, rarcurve,
                         composite_initial_point,
                         last_point_in_rarefaction,
                         odesolver,
                         deltaxi,
                         0, 0,
                         where_composite_begins, fam, 
                         new_rarcurve,
                         compositecurve,
                         final_direction,
                         reason_why,
                         edge);
        }

        virtual int curve(const AccumulationFunction *RarAccum, const FluxFunction *RarFlux,
                  const Boundary *RarBoundary, 
                  const Curve &rarcurve,
//                  std::vector<RealVector> &rarcurve, std::vector<double> &lambda,
                  const RealVector &composite_initial_point,
                  int last_point_in_rarefaction,
                  const ODE_Solver *odesolver,
                  double deltaxi,
                  void *linobj, double (*linear_function)(void *o, const RealVector &p),
                  int where_composite_begins, int fam, 
//                  std::vector<RealVector> &newrarcurve,
//                  std::vector<RealVector> &compositecurve,
                  Curve &new_rarcurve,
                  Curve &compositecurve,
                  RealVector &final_direction,
                  int &reason_why,
                  int &edge);

        virtual int correct_last_point(const ODE_Solver *odesolver, double deltaxi, WaveCurve &wavecurve);
                          
        static int double_contact_signal_event(const RealVector & where, double & determinant, int *obj, int * /*not used*/);
        static int rarefaction_of_composite_signal_event(const RealVector &where, double & current_diff_lambda, int *obj, int * /*not used*/);
        static int sigma_minus_lambda_signal_event(const RealVector &where, double &sigma_minus_lambda, int *obj, int * /*not used*/);

        // To be used only if use_maxsigma (above) is TRUE.
        //
        static int sigma_minus_maxsigma_signal_event(const RealVector &where, double &sigma_minus_maxsigma, int *obj, int * /*not used*/);

        static int explicit_bifurcation_expression_signal_event(const RealVector &where, double &expression, int *obj, int * /*not used*/);
        int transition_with_explicit_bifurcation(const ODE_Solver *odesolver, const RealVector &rarcmp_point, double init_time, RealVector &out, double &final_time);

//        static int characteristic_shock_signal_event(const RealVector &where, double &diff_lambda, int *obj, int * /*not used*/);

        
};

#endif // _COMPOSITECURVE_


#ifndef _RAREFACTIONCURVE_
#define _RAREFACTIONCURVE_

// Return values.
//
#define RAREFACTION_OK                           0
#define RAREFACTION_ERROR                        1
#define RAREFACTION_REACHED_INFLECTION           2
#define RAREFACTION_REACHED_BOUNDARY             3
#define RAREFACTION_INIT_OK                      4
#define RAREFACTION_INIT_ERROR                   5
#define RAREFACTION_COMPLEX_EIGENVALUE_AT_FAMILY 6
#define RAREFACTION_REACHED_COINCIDENCE_CURVE    7
#define RAREFACTION_REACHED_LINE                 8

// Options.
//
#define RAREFACTION_SPEED_SHOULD_INCREASE 10
#define RAREFACTION_SPEED_SHOULD_DECREASE 11

#define RAREFACTION    40
#define INTEGRAL_CURVE 41

#define RAREFACTION_INITIALIZE                   50
#define RAREFACTION_DONT_INITIALIZE              51

#include "ODE_Solver.h"
#include "Bisection.h"
#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Boundary.h"
#include "eigen.h"
#include "Utilities.h"

#include <complex> // To handle complex eigenpairs.



class WaveCurveFactory;

class RarefactionCurve {
    private:
    protected:
        const FluxFunction         *f;
        const AccumulationFunction *g;
        const Boundary             *b;

        const ODE_Solver           *s;

        int family;
        RealVector reference_vector;

        double directional_derivative(const RealVector &p, int fam, const RealVector &ref);

        void all_eigenvalues(const RealVector &p, int fam, std::vector<std::complex<double> > &lambda);
        void all_eigenvalues(const RealVector &p, int fam, RealVector &lambda);

        void add_point_to_curve(const RealVector &p, Curve &curve);

       
    public:
        RarefactionCurve(const AccumulationFunction *gg, const FluxFunction *ff, const Boundary *bb);
        virtual ~RarefactionCurve();

        int curve(const RealVector &initial_point,
                  int curve_family,
                  int increase,
                  int type_of_rarefaction, // For itself or as engine for integral curve.
                  int should_initialize,
                  const RealVector *direction,
                  const ODE_Solver *odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
                  double deltaxi,
                  Curve &rarcurve,
                  std::vector<RealVector> &inflection_points, // Will these survive/be added to the Curve class?
                  RealVector &final_direction,
                  int &reason_why, // Similar to Composite.
                  int &edge){

            return curve(initial_point,
                         curve_family,
                         increase,
                         type_of_rarefaction, // For itself or as engine for integral curve.
                         should_initialize,
                         direction,
                         odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
                         deltaxi,
                         0, 0, 
                         rarcurve,
                         inflection_points, // Will these survive/be added to the Curve class?
                         final_direction,
                         reason_why, // Similar to Composite.
                         edge);
        }

        int curve(const RealVector &initial_point,
                  int curve_family,
                  int increase,
                  int type_of_rarefaction, // For itself or as engine for integral curve.
                  int should_initialize,
                  const RealVector *direction,
                  const ODE_Solver *odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
                  double deltaxi,
                  void *linobj, double (*linear_function)(void *o, const RealVector &p),
                  Curve &rarcurve,
                  std::vector<RealVector> &inflection_points, // Will these survive/be added to the Curve class?
                  RealVector &final_direction,
                  int &reason_why, // Similar to Composite.
                  int &edge);

        int curve_from_boundary(const RealVector &initial_point, int side, 
                  int curve_family,
                  int increase,
                  int type_of_rarefaction, // For itself or as engine for integral curve.
                  const ODE_Solver *odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
                  double deltaxi,
                  Curve &rarcurve,
                  std::vector<RealVector> &inflection_points, // Will these survive/be added to the Curve class?
                  RealVector &final_direction,
                  int &reason_why, // Similar to Composite.
                  int &edge);

        static int field(int *neq, double *xi, double *in, double *out, int *obj, double* /* Not used */);
        static int Jacobian_field();

        static int inflection_signal_event(const RealVector & where, double & directional_derivative_measure, int *signal_object, int * /* reference_direction */);

        static int elliptic_region_signal_event_2D2D(const RealVector & where, double &discriminant, int *signal_object, int * /* reference_direction */);
        static int elliptic_region_signal_event_3D2D(const RealVector & where, double &discriminant, int *signal_object, int * /* reference_direction */);

        int initialize(const RealVector &p, int family, const RealVector &direction, RealVector &ref, double &dd);
        int initialize(const RealVector &p, int family, int increase, RealVector &ref, double &dd);

        friend class WaveCurveFactory;
        friend class ShockCurve;

       

//        static int curve(RarefactionCurve *obj, const RealVector &initial_point,
//                  int curve_family,
//                  int increase,
//                  int type_of_rarefaction, // For itself or as engine for integral curve.
//                  int should_initialize,
//                  const RealVector *direction,
//                  const ODE_Solver *odesolver, // Should it be another one for the Bisection? Can it really be const? If so, how to use initialize()?
//                  double deltaxi,
//                  std::vector<RealVector> &inflection_points, // Will these survive/be added to the Curve class?
//                  RealVector &final_direction,
//                  int &reason_why, // Similar to Composite.
//                  int &edge){

//            Curve rarcurve;
//            int info = obj->curve(intial_point, curve_family, increase, type_of_rarefaction, 
//                              should_initialize, direction, odesolver, 
//                              deltaxi, rarcurve, inflection_points, final_direction, reason_why, edge);

//            plot(rarcurve);
//            return info;
//        }
};

#endif // _RAREFACTIONCURVE_


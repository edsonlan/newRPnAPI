#ifndef _HUGONIOTCONTINUATION_
#define _HUGONIOTCONTINUATION_

#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "Boundary.h"
#include "RealVector.h"
#include "ReferencePoint.h"
#include "DoubleMatrix.h"
#include "eigen.h"
#include "JetMatrix.h"
#include "HugoniotCurve.h"

#ifndef HUGONIOTCONTINUATION_INITIALIZE_NO
#define HUGONIOTCONTINUATION_INITIALIZE_NO 0
#endif

#ifndef HUGONIOTCONTINUATION_INITIALIZE_YES
#define HUGONIOTCONTINUATION_INITIALIZE_YES 1
#endif

#ifndef HUGONIOTCONTINUATION_INITIALIZED_OK
#define HUGONIOTCONTINUATION_INITIALIZED_OK 10
#endif

#ifndef HUGONIOTCONTINUATION_INITIALIZE_ERROR
#define HUGONIOTCONTINUATION_INITIALIZE_ERROR 11
#endif

#ifndef HUGONIOTCONTINUATION_CURVE_OK
#define HUGONIOTCONTINUATION_CURVE_OK 20
#endif

#ifndef HUGONIOTCONTINUATION_CURVE_ERROR
#define HUGONIOTCONTINUATION_CURVE_ERROR 21
#endif

#ifndef HUGONIOTCONTINUATION_HYPERPLANE_OK
#define HUGONIOTCONTINUATION_HYPERPLANE_OK 30
#endif

#ifndef HUGONIOTCONTINUATION_HYPERPLANE_ERROR
#define HUGONIOTCONTINUATION_HYPERPLANE_ERROR 31
#endif

#ifndef HUGONIOTCONTINUATION_DIRECTION_OK
#define HUGONIOTCONTINUATION_DIRECTION_OK 40
#endif

#ifndef HUGONIOTCONTINUATION_DIRECTION_ERROR
#define HUGONIOTCONTINUATION_DIRECTION_ERROR 41
#endif

#ifndef HUGONIOTCONTINUATION_NEWTON_OK
#define HUGONIOTCONTINUATION_NEWTON_OK 50
#endif

#ifndef HUGONIOTCONTINUATION_NEWTON_ERROR
#define HUGONIOTCONTINUATION_NEWTON_ERROR 51
#endif

#ifndef HUGONIOTCONTINUATION_NEWTON_STEP_OK
#define HUGONIOTCONTINUATION_NEWTON_STEP_OK 60
#endif

#ifndef HUGONIOTCONTINUATION_NEWTON_STEP_ERROR
#define HUGONIOTCONTINUATION_NEWTON_STEP_ERROR 61
#endif

#ifndef HUGONIOTCONTINUATION_SET_BIFURCATION_PLANE_OK
#define HUGONIOTCONTINUATION_SET_BIFURCATION_PLANE_OK 70
#endif

#ifndef HUGONIOTCONTINUATION_SET_BIFURCATION_PLANE_ERROR
#define HUGONIOTCONTINUATION_SET_BIFURCATION_PLANE_ERROR 71
#endif

#ifndef HUGONIOTCONTINUATION_NEAR_COINCIDENCE_CURVE
#define HUGONIOTCONTINUATION_NEAR_COINCIDENCE_CURVE 81
#endif

#ifndef SQRT_TWO
#define SQRT_TWO 1.41421356237
#endif

#ifndef MAX_COS_ANGLE
#define MAX_COS_ANGLE (0.96)
#endif

// This class computes the Hugoniot locus using a continuation method.
// 
// The class is abstract, as the following methods must be supplied by derived classes:
//
//     void jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH); 
//     int fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &hyperplane, RealVector &Hugoniot_direction);
//

// TODO. If the necessity ever arises of having a copy constructor or assignment operator, refer to
//
//     http://www.parashift.com/c++-faq-lite/copy-of-abc-via-clone.html
//
class HugoniotContinuation: public HugoniotCurve {
    private:
    protected:
        const FluxFunction *f;
        const AccumulationFunction *g;
        const Boundary *b;

        ReferencePoint ref;
        
        double bifurcation_space_coordinate;
        int bifurcation_space_coordinate_index;
        bool there_is_a_bifurcation_space;

        // 
        double default_step_size_;
        double max_default_step_size_;

        // Derived classes may modify the following:
        //
        virtual int find_initial_direction(const RealVector &in, const RealVector &hint_direction, int fam, RealVector &initial_direction);

        virtual int fill_hyperplane(const RealVector &origin, DoubleMatrix &hyperplane);

        virtual void Newton_system_for_Hugoniot(const RealVector &p, const DoubleMatrix &hyperplane, DoubleMatrix &Newton_matrix, RealVector &error);

        virtual int Newton_in_hyperplane(const RealVector &origin, const DoubleMatrix &hyperplane, RealVector &Hugoniot_intersection);

        // TODO: Decide which version of the following methods is to survive. Right now (16-Aug-2013) the version that doesn't write on hyperplane
        //       works fine.
        virtual int Newton_step(const RealVector &previous_point, double &step_size, int &number_of_steps_with_unchanged_size, const RealVector &previous_direction,
                                DoubleMatrix &hyperplane, RealVector &Hugoniot_intersection);

        virtual int Newton_step(const RealVector &previous_point, double &step_size, int &number_of_steps_with_unchanged_size, const RealVector &previous_direction,
                                RealVector &Hugoniot_intersection);

        virtual int fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &hyperplane, RealVector &Hugoniot_direction);
        
        // Derived classes must provide:
        //
        virtual RealVector orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane) = 0;
    public:
        HugoniotContinuation(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb);
        virtual ~HugoniotContinuation();

        // Derived classes must provide:
        //
        // virtual void jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH) = 0; // Is no more.

        // F, JF, G and JG must be computed before entering here.
        // TODO: Pass a flag to tell if it must compute H or H and nablaH.
        //
        virtual void jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                                  const RealVector &G, const DoubleMatrix &JG, 
                                  RealVector &H, DoubleMatrix &nablaH){
            RealVector C(0);
            DoubleMatrix JC(0, 0);

            jet_Hugoniot(F, JF, G, JG, C, JC, H, nablaH);

            return;
        }

        virtual void jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                                  const RealVector &G, const DoubleMatrix &JG, 
                                  const RealVector &C, const DoubleMatrix &JC, 
                                  RealVector &H, DoubleMatrix &nablaH) = 0;  

        virtual void set_reference_point(const ReferencePoint &r);

        virtual int set_bifurcation_space_coordinate(int Theta_index);
        virtual bool reference_point_near_coincidence();

        // Needed by ShockCurve.
        //
        virtual int find_continuation_direction(const RealVector &Hugoniot_point, const RealVector &hint, RealVector &Hugoniot_direction);

        // Shockspeed between two generic points.
        // 
        virtual double shockspeed(const FluxFunction *fp, const AccumulationFunction *gp, const RealVector &Up,
                                  const FluxFunction *fm, const AccumulationFunction *gm, const RealVector &Um) const;

        // Shockspeed (and its derivative) between point U and the reference point. 
        //
        // TODO: A generic version could also be written.
        //
        virtual void jet_sigma(const RealVector &U, const RealVector &Hugoniot_direction, double &sigma, double &sigma_dot) const;

        // Hugoniot_direction is only needed when computing the derivative.
        // Therefore, when computing only sigma it can be passed a null pointer instead.
        //
        virtual void sigma_jet(const RealVector &U, const RealVector *Hugoniot_direction, int degree, JetMatrix &s) const;

        virtual double sigma(const RealVector &F, const RealVector &G) const;

        virtual double sigma(const RealVector &Fp, const RealVector &Gp, const RealVector &Fm, const RealVector &Gm) const;

        virtual int curve_point(const RealVector &previous_point, double previous_sigma_between_points,
                                const RealVector &direction, 
                                int &step_size_increased,
                                double &step_size, int &number_of_steps_with_unchanged_size, 
                                RealVector &Hugoniot_intersection, double &sigma_between_points,
                                RealVector &Hugoniot_direction);

        virtual int curve_engine(const RealVector &in, const RealVector &initial_direction, 
                                 RealVector &final_direction, std::vector<RealVector> &shockcurve, int &edge);

        virtual int curve(std::vector< std::vector<RealVector> > &curve);
        virtual int disconnected_curve(const RealVector &point_on_edge, int entry_edge, std::vector<RealVector> &curve, int &edge);
        virtual int all_curves(std::vector< std::vector<RealVector> > &curves, std::vector<int> &edges);
        virtual int connected_curve(std::vector< std::vector<RealVector> > &curves, std::vector<int> &edges);

        virtual const FluxFunction         * flux()         const {return f;}
        virtual const AccumulationFunction * accumulation() const {return g;}
        virtual const Boundary             * boundary()     const {return b;}

        // Get/set default_step_size and max_default_step_size
        virtual double default_step_size(){return default_step_size_;}

        virtual void default_step_size(double d){
            default_step_size_ = d;
            max_default_step_size_ = 100.0*default_step_size_;
            return;
        }

        // TODO: Move the code below to a Stone-something ASAP.
        void hugoniot_on_side(double p, int side, double &RHv, double &RHv_prime);
        int Newton_on_side(double init_p, int side, double &p);

        bool find_a_point_on_a_side(int side, RealVector &p);
        bool find_a_point_on_a_side(int side, 
                                                  const std::vector<RealVector> &side_points,
                                                  const std::vector<bool> &use_this_segment,
                                                  RealVector &p);

        // These are needed to complete the inheritance from HugoniotCurve.
        //
        virtual void list_of_reference_points(std::vector<int> &type, std::vector<std::string> &name) const {
            name.clear();
            name.push_back(std::string("Generic point"));

            type.clear();
            type.push_back(0);

            return;
        } 

        virtual void curve(const ReferencePoint &ref, int type, std::vector<Curve> &c){
            set_reference_point(ref);

            std::vector< std::vector<RealVector> > curves;
            std::vector<int> edges;

            connected_curve(curves, edges);

            c.clear();
            for (int i = 0; i < curves.size(); i++){
                Curve ctemp;
                ctemp.curve = curves[i];

                c.push_back(ctemp);
            }

            return;
        }
};

#endif // _HUGONIOTCONTINUATION_


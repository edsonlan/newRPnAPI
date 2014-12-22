#ifndef EXPLICIT_BIFURCATION_CURVES
#define EXPLICIT_BIFURCATION_CURVES

#include <vector>
#include "RealVector.h"
#include "Curve.h"

#define EXPLICIT_BIFURCATION_CODE_OK             0
#define EXPLICIT_BIFURCATION_CODE_ERROR          1

#define SHOCK_CROSSED_EXPLICIT_BIFURCATION       10
#define SHOCK_DID_NOT_CROSS_EXPLICIT_BIFURCATION 11

class Explicit_Bifurcation_Curves {
    private:
    protected:
    public:
        Explicit_Bifurcation_Curves();
        virtual ~Explicit_Bifurcation_Curves();

        // Explicit secondary bifurcation curve.
        //
        virtual void expl_sec_bif_crv(int vertex, int nos, 
                                      std::vector<RealVector> &vertex_to_umbilic, 
                                      std::vector<RealVector> &umbilic_to_side) = 0;

//        // The shock curve crossed the secondary bifurcation.
//        //
//        virtual int cross_sec_bif(const RealVector &previous_point, const RealVector &point, RealVector &crossing_point, int &region) = 0;

        // Find the correspondences on the secondary bifurcations
        //
        virtual int sec_bif_correspondence(int side_opposite_vertex, int nos, 
                                           std::vector<RealVector> &point, 
                                           std::vector<RealVector> &correspondent_point,
                                           int &model_specific_error_code) = 0;

        // Compute the expressions that form the transition lines.
        //
        virtual RealVector expressions(const RealVector &point) = 0;

        virtual void subdivide_curve_in_regions(const Curve &curve, 
                                        Curve &classified_curve,
                                        std::vector<int> &Lax_transition_index,
                                        std::vector<int> &region_transition_index) = 0;

        virtual void subdivide_segmented_curve_in_regions(const std::vector<RealVector>&curve, 
                                                  std::vector<RealVector> &classified_curve,
                                                  std::vector<int> &Lax_transition_index,
                                                  std::vector<int> &region_transition_index) = 0;
};

#endif // EXPLICIT_BIFURCATION_CURVES


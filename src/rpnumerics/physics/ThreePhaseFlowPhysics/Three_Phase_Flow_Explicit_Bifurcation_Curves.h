#ifndef _THREE_PHASE_FLOW_EXPLICIT_BIFURCATION_CURVES_
#define _THREE_PHASE_FLOW_EXPLICIT_BIFURCATION_CURVES_

#include "Explicit_Bifurcation_Curves.h"

#include "Curve.h"
//#include "Three_Phase_FluxFunction.h"
#include "ThreePhaseFlowFluxFunction.h"
#include "IsoTriang2DBoundary.h"
#include "Utilities.h"

#define PERMEABILITY_PARAMETERS_OK    0
#define PERMEABILITY_PARAMETERS_ERROR 1

#define FLUX_PARAMETERS_OK    2
#define FLUX_PARAMETERS_ERROR 3
 
// The domain is divided in regions according
// to the value of the inequalities that define the explicit bifurcations curves.
//
// Let the triangle below be divided in the following regions:
//
//     The line between the "pure W" vertex and the "No W" side divides the triangle in two regions:
//     below (W-) and above (W+) it. 
//
//     The line between the "pure O" vertex and the "No O" side divides the triangle in two regions:
//     to the left (O-) and to the right (O+) of it.
//
//     The line between the "pure G" vertex and the "No G" side divides the triangle in two regions:
//     below (G-) and above (G+) it. 
//
// Any point in the domain can be classified according to where in all of the three abovementioned regions is.
// The suffixes "-" and "+" above are replaced by "M" and "P" respectively.
//
//
//             pure O
//               /|\
//              / | \
//             /  |  \
//            /   |   \
//     No W  /    |    \ No G
//          /     |     \
//         /  -   |  +   \
//        /       |       \
//       /        |        \
// pure G-------------------pure W
//              No O
// The 

#define REGION_WM_OM_GM 0
#define REGION_WM_OP_GM 1
#define REGION_WP_OP_GM 2
#define REGION_WP_OP_GP 3
#define REGION_WP_OM_GP 4
#define REGION_WM_OM_GP 5

class Three_Phase_Flow_Explicit_Bifurcation_Curves : public Explicit_Bifurcation_Curves {
    private:
    protected:
//        Three_Phase_FluxFunction *flux;
        ThreePhaseFlowFluxFunction *flux;

//        void points_to_segments(std::vector<RealVector> &v);

        virtual RealVector sec_bif_correspondence(int side_opposite_vertex, const RealVector &point, const RealVector &mu); /* Ok */
    public:
        Three_Phase_Flow_Explicit_Bifurcation_Curves(ThreePhaseFlowFluxFunction *f); /* Ok */
        virtual ~Three_Phase_Flow_Explicit_Bifurcation_Curves(); /* Ok */

        virtual void vertex_and_side(int side_opposite_vertex, const RealVector &mu, RealVector &vertex, RealVector &point_on_side); /* Ok */

        // Concrete method declared in the parent class Explicit_Bifurcation_Curves.
        // 
        int sec_bif_correspondence(int side_opposite_vertex, int nos, 
                                   std::vector<RealVector> &point, 
                                   std::vector<RealVector> &correspondent_point,
                                   int &model_specific_error_code); /* Ok */

        // Concrete method declared in the parent class Explicit_Bifurcation_Curves.
        //
        RealVector expressions(const RealVector &point); /* Ok */

        int region(const RealVector &expressions); /* Ok */

        void subdivide_curve_in_regions(const Curve &curve, 
                                        Curve &classified_curve,
                                        std::vector<int> &Lax_transition_index,
                                        std::vector<int> &region_transition_index); // TODO: Add an "int region" member to Curve. /* Ok */

        void subdivide_segmented_curve_in_regions(const std::vector<RealVector>&curve, 
                                                  std::vector<RealVector> &classified_curve,
                                                  std::vector<int> &Lax_transition_index,
                                                  std::vector<int> &region_transition_index);

        /* Ok */
        virtual int check_permeability_parameters(int &model_specific_error_code);

        /* Ok */
        virtual int check_flux_parameters(int &model_specific_error_code);
};

#endif // _THREE_PHASE_FLOW_EXPLICIT_BIFURCATION_CURVES_


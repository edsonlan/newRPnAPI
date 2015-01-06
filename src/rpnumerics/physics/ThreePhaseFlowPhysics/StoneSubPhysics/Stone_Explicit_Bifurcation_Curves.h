#ifndef STONE_EXPLICIT_BIFURCATION_CURVES
#define STONE_EXPLICIT_BIFURCATION_CURVES

#include "Explicit_Bifurcation_Curves.h"
#include "StoneFluxFunction.h"
#include "Three_Phase_Boundary.h"
#include "Three_Phase_Flow_Explicit_Bifurcation_Curves.h"

class Stone_Explicit_Bifurcation_Curves : public Three_Phase_Flow_Explicit_Bifurcation_Curves {
    private:
    protected:
        StoneFluxFunction *f;

        void line(const RealVector &p, const RealVector &q, int nos, std::vector<RealVector> &v); /* OK */
//        RealVector sec_bif_correspondence(int side_opposite_vertex, const RealVector &point, const RealVector &mu);

//        void vertex_and_side(int side_opposite_vertex, const RealVector &mu, RealVector &vertex, RealVector &point_on_side);

    public:
        Stone_Explicit_Bifurcation_Curves(StoneFluxFunction *ff); /* OK */
        virtual ~Stone_Explicit_Bifurcation_Curves(); /* OK */

        /* OK */
        // Concrete method declared in the parent class Explicit_Bifurcation_Curves.
        //
        void expl_sec_bif_crv(int side_opposite_vertex, int nos, 
                              std::vector<RealVector> &vertex_to_umbilic, 
                              std::vector<RealVector> &umbilic_to_side);

//        RealVector equations(const RealVector &point);
//        int region(const RealVector &equations);

//        virtual int cross_sec_bif(const RealVector &previous_point, const RealVector &point, RealVector &crossing_point, int &region);

//        int sec_bif_correspondence(int side_opposite_vertex, int nos, 
//                                   std::vector<RealVector> &point, 
//                                   std::vector<RealVector> &correspondent_point,
//                                   int &model_specific_error_code);

        int check_permeability_parameters(int &model_specific_error_code); /* OK */
        int check_flux_parameters(int &model_specific_error_code); /* OK */
};

#endif // STONE_EXPLICIT_BIFURCATION_CURVES


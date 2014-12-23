#include "Boundary.h"

#include "Extension_Curve.h"
#include "Envelope_Curve.h"

#include "VLE_Flash_TPCW.h"

#define DOMAIN_IS_VAPOR  0
#define DOMAIN_IS_LIQUID 1

#define DOMAIN_EDGE_LOWER_COMPOSITION 10
#define DOMAIN_EDGE_UPPER_COMPOSITION 20

// TODO: Explain what these names mean.
#define SINGLEPHASE_TOTAL_COMPOSITION_SIDE 30
#define SINGLEPHASE_ANY_OTHER_SIDE         40

// Sides
#define SINGLEPHASE_LOWER_COMPOSITION   100
#define SINGLEPHASE_MINIMUM_TEMPERATURE 200
#define SINGLEPHASE_UPPER_COMPOSITION   300
#define SINGLEPHASE_MAXIMUM_TEMPERATURE 400

class SinglePhaseBoundary : public Boundary {
    private:
    protected:
        RealVector bmin, bmax;

        // This is the value of the composition that is above zero and below the maximum value of the composition.
        // This value will be used by inside().
        //
        double top_min_composition; 

        VLE_Flash_TPCW *Flash;

        int domain_type;

        void init(VLE_Flash_TPCW *f, double Theta_min, double Theta_max, int domain_type_, double (*d)(double Theta));
        void copy_from(const SinglePhaseBoundary &orig);

        int edge_segments(int where_constant, int number_of_steps, std::vector<RealVector> &seg);

        static double default_dimensionalize_variable(double Theta);
        double (*dimensionalize_variable)(double Theta);

        int top_boundary_points;
    public:
        SinglePhaseBoundary(VLE_Flash_TPCW *f, double Theta_min, double Theta_max, int domain_type_);
        SinglePhaseBoundary(VLE_Flash_TPCW *f, double Theta_min, double Theta_max, int domain_type_, double (*d)(double Theta));
        SinglePhaseBoundary(const SinglePhaseBoundary &orig);
        SinglePhaseBoundary(const SinglePhaseBoundary *orig);

        ~SinglePhaseBoundary();

        bool inside(const RealVector &y) const;
        bool inside(const double *) const;

        //! Virtual constructor
        Boundary * clone()const;

        //! Minimums boundary values accessor
        const RealVector & minimums() const;

        //! Maximums boundary values accessor
        const RealVector & maximums() const;

        RealVector intersect(RealVector &y1, RealVector &y2) const;

        //! Returns the boundary type
        const char * boundaryType() const;

        int intersection(const RealVector &p, const RealVector &q, RealVector &r, int &) const;

        void extension_curve(const FluxFunction *f, const AccumulationFunction *a,
                             GridValues &gv,
                             int where_constant, int number_of_steps, bool singular,
                             int fam, int characteristic,
                             std::vector<RealVector> &c, std::vector<RealVector> &d);

        void extension_curve(const FluxFunction *df, const AccumulationFunction *da, // Over the domain
                             const FluxFunction *cf, const AccumulationFunction *ca, // Over the curve 
                             GridValues &gv,
                             int where_constant, int number_of_steps, bool singular,
                             int fam, int characteristic,
                             std::vector<RealVector> &c, std::vector<RealVector> &d);

        void envelope_curve(const FluxFunction *f, const AccumulationFunction *a,
                            GridValues &gv,
                            int where_constant, int number_of_steps, bool singular,
                            std::vector<RealVector> &c, std::vector<RealVector> &d);

        void set_number_of_top_boundary_points(int p){top_boundary_points = p; return;}

        RealVector side_transverse_interior(const RealVector &p, int s) const;
        RealVector side_tangent_vector(const RealVector &p, int s) const;

        //void set_dimensionalization(double (*f)(double Theta)){dimensionalize_variable = f; return;}

        void edge_segments(int where_constant, int number_of_steps, std::vector<RealVector> &seg) const;
        void list_of_sides(std::vector<int> &where_constant_codes, std::vector<std::string> &where_constant_names) const;
        void physical_boundary(std::vector<std::vector<RealVector> > &pb) const;
};


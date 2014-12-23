#ifndef _COLORCURVE_
#define _COLORCURVE_

#include <vector>
#include <deque>
#include <string>

#include "FluxFunction.h"
#include "AccumulationFunction.h"

#include "eigen.h"
#include "ViscosityJetMatrix.h"
#include "Viscosity_Matrix.h"
#include "ReferencePoint.h"

#define INTERPOLATION_ERROR -1
#define INTERPOLATION_OK     0

//#define _SHOCK_SIMPLE_ACCUMULATION_  10  // Traditional rarefaction, using dgeev.
//#define _SHOCK_GENERAL_ACCUMULATION_ 11  // Rarefaction with generalized eigenpairs, using dggev.

#ifndef UNCLASSIFIABLE_POINT
#define UNCLASSIFIABLE_POINT (-1)
#endif

extern "C" void dgesv_(int*, int*, double *, int*, int *, double *, int*, int*);

struct HugoniotPolyLine {
public:
    // Elements of extrema, points coordinates, speeds and eginvalues
    //
    std::vector<RealVector> point; // Each element has size = dimension (dim).
    std::vector<double> speed; // Speed at each point.
    std::vector<RealVector> eigenvalue; // Each element has size = number of valid eigenvalues (fam).
    // Notice that [ mu ~= lambda - s ] is stored.

    // Elements of segment, type (color) and signature
    //
    std::vector<int> type; // A color table is needed for graphical issues.
    std::vector<std::string> signature; // This returns the Hugoniot signature, i.e., "--++", etc
    // zero means characteristic, "-." or "+." means complex
    // conjugate, with the real part sign.

    HugoniotPolyLine() {
    };

    ~HugoniotPolyLine() {
    };
};

// Color table. Could be anywhere.
//struct ColorTable{
//    public:
//        static int color[16][3];
//};

class ColorCurve {
private:
    std::string sp, sm, sc, sz;
protected:
    int solve(const double *A, double *b, int dim, double *x);

    void Left_Newton_improvement(const RealVector &input, const int type, RealVector &out);
    void Right_Newton_improvement(const RealVector &input, const int type, RealVector &out);

    int interpolate(const RealVector &p, double &s_p,
            const std::vector<double> &eigenvalue_p, const int type_p,
            const RealVector &q, double &s_q,
            const std::vector<double> &eigenvalue_q, const int type_q,
            vector<RealVector> &r, vector<int> &rtype);

    int complete_point(RealVector &p, double &s, std::vector<double> &eigenvalue, int *complex);

    int classify_point(RealVector &p, double &s, std::vector<double> &eigenvalue, std::string &signature);

    void classify_segment(RealVector &p, RealVector &q,
            std::vector<HugoniotPolyLine> &classified_curve,
            std::vector<RealVector> &transition_list);

    void classify_segment_with_data(
            RealVector &p, double &s_p, std::vector<double> &eigenvalue_p, std::string &ct_p, int &type_p,
            RealVector &q, double &s_q, std::vector<double> &eigenvalue_q, std::string &ct_q, int &type_q,
            HugoniotPolyLine &hpl,
            std::vector<HugoniotPolyLine> &classified_curve,
            std::vector<RealVector> &transition_list);

    const FluxFunction * fluxFunction_;
    const AccumulationFunction * accFunction_;

    const Viscosity_Matrix *vm;

    RealVector ref_point;
    std::vector<double> ref_eigenvalue;
    std::vector<double> ref_e_complex;

    // The values at reference point are stored, as the Viscosity Matrix
    RealVector F_ref, G_ref;
    DoubleMatrix JF_ref, JG_ref;
    ViscosityJetMatrix B_ref;
public:

    ColorCurve(const FluxFunction &, const AccumulationFunction &, const Viscosity_Matrix *v = 0);
//    ColorCurve(const FluxFunction *, const AccumulationFunction *);
    virtual ~ColorCurve();

    void classify_segmented_curve(std::vector<RealVector> &original, const ReferencePoint &ref,
            std::vector<HugoniotPolyLine> &classified_curve,
            std::vector<RealVector> &transition_list);

    void classify_continuous_curve(std::deque<RealVector> &original, const ReferencePoint &ref,
            HugoniotPolyLine &classified_curve,
            std::vector<RealVector> &transition_list);
};

#endif // _COLORCURVE_

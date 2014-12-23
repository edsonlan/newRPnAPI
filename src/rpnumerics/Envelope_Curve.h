#ifndef _ENVELOPE_CURVE_
#define _ENVELOPE_CURVE_

#include "TwoImplicitFunctions.h"
#include "Contour2p5_Method.h"
#include "GridValues.h"

class Envelope_Curve : public TwoImplicitFunctions {
private:
protected:
    Matrix<double> segment_flux, segment_accum, segment_data;
        
    const FluxFunction *ff;
    const AccumulationFunction *aa;
public:
    // TODO: Maybe this class could be formed by purely static. In that case the ctor() may be useless.
    // The convenience of this approach is to be discussed sometime.
    Envelope_Curve(){
        gv = 0;
        oc = 0;
        singular = true;
        
        segment_flux.resize(2, 2);
        segment_accum.resize(2, 2);
        
        segment_data.resize(4, 2);
    }
    
    ~Envelope_Curve(){}

    bool valid_segment(int i);

    int function_on_vertices(double *foncub, int domain_i, int domain_j, int kl);
        
    void curve(const FluxFunction *f, const AccumulationFunction *a, 
               GridValues &g, bool is_singular, 
               std::vector<RealVector> &original_curve,
               std::vector<RealVector> &envelope_on_curve, std::vector<RealVector> &envelope_on_domain);
};

#endif // _ENVELOPE_CURVE_


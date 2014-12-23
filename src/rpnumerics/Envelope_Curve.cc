#include "Envelope_Curve.h"

int Envelope_Curve::function_on_vertices(double *foncub, int domain_i, int domain_j, int kl) {
    if (gv == 0 || oc == 0) return INVALID_FUNCTION_ON_VERTICES;

    // Some differences of the Rankine-Hugoniot relation
    double dg0 = segment_accum(0, kl) - gv->G_on_grid(domain_i, domain_j).component(0);
    double dg1 = segment_accum(1, kl) - gv->G_on_grid(domain_i, domain_j).component(1);
    double df0 = segment_flux(0, kl)  - gv->F_on_grid(domain_i, domain_j).component(0);
    double df1 = segment_flux(1, kl)  - gv->F_on_grid(domain_i, domain_j).component(1);

    // RH relation, without shock speed
    foncub[0] = dg1 * df0 - dg0 * df1;

    // dRH/ds, with s the "arc length" of the segment.
    foncub[1] = dg1 * segment_data(0 + kl) - dg0 * segment_data(2 + kl)
              + segment_data(6 + kl) * df0 - segment_data(4 + kl) * df1;

    return VALID_FUNCTION_ON_VERTICES;
}

bool Envelope_Curve::valid_segment(int i){
    if (oc == 0) return false;

    double JF[2][2], JG[2][2], F[2], G[2];

    double d_eta[2];
    for (int j = 0; j < 2; j++) d_eta[j] = oc->at(i).component(j) - oc->at(i + 1).component(j);
    
    for (int j = 0; j < 2; j++){
        ff->fill_with_jet(2, oc->at(i + j).components(), 1, F, &JF[0][0], 0);
        aa->fill_with_jet(2, oc->at(i + j).components(), 1, G, &JG[0][0], 0);
       
        segment_data(0 + j) = JF[0][0]*d_eta[0] + JF[0][1]*d_eta[1];
        segment_data(2 + j) = JF[1][0]*d_eta[0] + JF[1][1]*d_eta[1];
        segment_data(4 + j) = JG[0][0]*d_eta[0] + JG[0][1]*d_eta[1];
        segment_data(6 + j) = JG[1][0]*d_eta[0] + JG[1][1]*d_eta[1];

        for (int k = 0; k < 2; k++){
            segment_flux(k, j)  = F[k];
            segment_accum(k, j) = G[k];
        }
    }

    return true;
}

void Envelope_Curve::curve(const FluxFunction *f, const AccumulationFunction *a, 
                           GridValues &g, bool is_singular, 
                           std::vector<RealVector> &original_curve,
                           std::vector<RealVector> &envelope_on_curve,
                           std::vector<RealVector> &envelope_on_domain){
                            
    ff = f;
    aa = a;
    
    singular = is_singular;

    gv = &g;
    oc = &original_curve;
    
    gv->fill_functions_on_grid(ff, aa);
    
    envelope_on_curve.clear(); 
    envelope_on_domain.clear(); 

    Contour2p5_Method::contour2p5(this, envelope_on_curve, envelope_on_domain);
    
    return;
}


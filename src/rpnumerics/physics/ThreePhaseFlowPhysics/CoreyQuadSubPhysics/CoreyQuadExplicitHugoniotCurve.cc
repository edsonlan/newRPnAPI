#include "CoreyQuadExplicitHugoniotCurve.h"
#include "CoreyQuadSubPhysics.h"

CoreyQuadExplicitHugoniotCurve::CoreyQuadExplicitHugoniotCurve(CoreyQuadSubPhysics *c) : HugoniotCurve(c->flux(), c->accumulation(), c->boundary()), corey(c) {
    side_opposite_vertex.push_back(THREE_PHASE_BOUNDARY_SW_ZERO);
    side_opposite_vertex.push_back(THREE_PHASE_BOUNDARY_SO_ZERO);
    side_opposite_vertex.push_back(THREE_PHASE_BOUNDARY_SG_ZERO);

    method_ = EXPLICIT_HUGONIOT;

    info_ = std::string("CoreyQuad Explicit");
}

CoreyQuadExplicitHugoniotCurve::~CoreyQuadExplicitHugoniotCurve(){
}

// Compute the umbilic point given that grw == grg && grg == gro.
//
RealVector CoreyQuadExplicitHugoniotCurve::compute_umbilic_point(const RealVector &mu){
    double inv_sum = 1.0/sum(mu);
    RealVector umbilic_point(2);
    for (int i = 0; i < 2; i++) umbilic_point(i) = mu(i)*inv_sum;

    return umbilic_point;
}

//void CoreyQuadExplicitHugoniotCurve::curve(const ReferencePoint &ref, int type, std::vector<HugoniotPolyLine> &c){
void CoreyQuadExplicitHugoniotCurve::curve(const ReferencePoint &ref, int type, std::vector<Curve> &c){
    c.clear();

//    RealVector flux_params = flux->fluxParams().params();

    double grw = corey->grw()->value();
    double grg = corey->grg()->value();
    double gro = corey->gro()->value();

    // Only works in this case:
    //
    if (grw == grg && grg == gro){
        // Parameters that will be used.
        //
        muw = corey->muw()->value();
        mug = corey->mug()->value();
        muo = corey->muo()->value();

        RealVector mu(3);
        mu(0) = muw; // muw
        mu(2) = mug; // mug
        mu(1) = muo; // muo

        // Project the reference point according to the case.
        //
//        reference_point.point = project(ref.point, type);
        reference_point = ref;

        int n = 100;

        switch (type){
            case COREYQUAD_GENERIC_POINT:
            {
                double phi_begin = 0.0;
                double phi_end   = M_PI;

                std::vector<Curve> curve;
                ParametricPlot::plot(&generic, &f_asymptote, (void*)this, phi_begin, phi_end, n, boundary, curve);
                for (int i = 0; i < curve.size(); i++) c.push_back(curve[i]);

                break;
            }
            case COREYQUAD_G_VERTEX:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }                

                // Add GW side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), point_on_side[1], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[1], corey->W(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                // Add GO side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), point_on_side[0], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[0], corey->O(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                // Add G-point_on_WO-side as two Curves.
                //
                {
                    RealVector umbilic_point = compute_umbilic_point(mu);

                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[2], n, temp_curve);
                    c.push_back(temp_curve);
                }

                break;
            }
            case COREYQUAD_W_VERTEX:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }                

                // Add GW side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), point_on_side[1], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[1], corey->W(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                // Add WO side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), point_on_side[2], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[2], corey->O(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                // Add W-point_on_WO-side as two Curves.
                //
                {
                    RealVector umbilic_point = compute_umbilic_point(mu);

                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[0], n, temp_curve);
                    c.push_back(temp_curve);
                }

                break;
            }
            case COREYQUAD_O_VERTEX:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }                

                // Add GO side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), point_on_side[0], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[0], corey->O(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                // Add WO side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), point_on_side[2], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[2], corey->O(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                // Add W-point_on_WO-side as two Curves.
                //
                {
                    RealVector umbilic_point = compute_umbilic_point(mu);

                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->O(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[1], n, temp_curve);
                    c.push_back(temp_curve);
                }

                break;
            }
            case COREYQUAD_GW_SIDE:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }                

                // Add GW side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), point_on_side[1], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[1], corey->W(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                double phi_begin = -M_PI/2 + .1;
                double phi_end   =  M_PI/2 + .1;

                std::vector<Curve> curve;
                ParametricPlot::plot(&GW_side, &f_asymptote, (void*)this, phi_begin, phi_end, n, boundary, curve);
                for (int i = 0; i < curve.size(); i++) c.push_back(curve[i]);

                break;
            }
            case COREYQUAD_WO_SIDE:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }                

                // Add WO side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), point_on_side[2], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[2], corey->O(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                double phi_begin = -M_PI/2;
                double phi_end   =  M_PI/2;

                std::vector<Curve> curve;
                ParametricPlot::plot(&WO_side, &f_asymptote, (void*)this, phi_begin, phi_end, n, boundary, curve);
                for (int i = 0; i < curve.size(); i++) c.push_back(curve[i]);

                break;
            }
            case COREYQUAD_GO_SIDE:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }                

                // Add GO side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), point_on_side[0], n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(point_on_side[0], corey->O(), n, temp_curve);
                    c.push_back(temp_curve);
                }

                double phi_begin = -M_PI/2;
                double phi_end   =  M_PI/2;

                std::vector<Curve> curve;
                ParametricPlot::plot(&GO_side, &f_asymptote, (void*)this, phi_begin, phi_end, n, boundary, curve);
                for (int i = 0; i < curve.size(); i++) c.push_back(curve[i]);

                break;
            }
            case COREYQUAD_G_BIFURCATION:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }

                RealVector umbilic_point = compute_umbilic_point(mu);

                // Add G-point_on_side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[2], n, temp_curve);
                    c.push_back(temp_curve);
                }            

                double phi_begin = 0.0;
                double phi_end   = M_PI;

                std::vector<Curve> curve;
                ParametricPlot::plot(&G_bif, &f_asymptote, (void*)this, phi_begin, phi_end, n, boundary, curve);
                for (int i = 0; i < curve.size(); i++) c.push_back(curve[i]);

                break;
            }
            case COREYQUAD_W_BIFURCATION:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }

                RealVector umbilic_point = compute_umbilic_point(mu);

                // Add W-point_on_side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[0], n, temp_curve);
                    c.push_back(temp_curve);
                }            

                double phi_begin = 0.0;
                double phi_end   = M_PI;

                std::vector<Curve> curve;
                ParametricPlot::plot(&W_bif, &f_asymptote, (void*)this, phi_begin, phi_end, n, boundary, curve);
                for (int i = 0; i < curve.size(); i++) c.push_back(curve[i]);

                break;
            }
            case COREYQUAD_O_BIFURCATION:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }

                RealVector umbilic_point = compute_umbilic_point(mu);

                // Add O-point_on_side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->O(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[1], n, temp_curve);
                    c.push_back(temp_curve);
                }            

                double phi_begin = 0.0;
                double phi_end   = M_PI;

                std::vector<Curve> curve;
                ParametricPlot::plot(&O_bif, &f_asymptote, (void*)this, phi_begin, phi_end, n, boundary, curve);
                for (int i = 0; i < curve.size(); i++) c.push_back(curve[i]);

                break;
            }
            case COREYQUAD_UMBILIC_POINT:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }

                RealVector umbilic_point = compute_umbilic_point(mu);

                // Add W-point_on_side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[0], n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add O-point_on_side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->O(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[1], n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add G-point_on_side as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, point_on_side[2], n, temp_curve);
                    c.push_back(temp_curve);
                }  

                break;
            }
            case COREYQUAD_E_POINT:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }

                RealVector umbilic_point = compute_umbilic_point(mu);
                RealVector E = point_on_side[0];
                RealVector B = point_on_side[1];
                RealVector D = point_on_side[2];

                // Add W-E as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, E, n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add O-E-G as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->O(), E, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(E, corey->G(), n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add B-D as one Curve.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(B, D, n, temp_curve);
                    c.push_back(temp_curve);
                }  

                break;
            }
            case COREYQUAD_B_POINT:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }

                RealVector umbilic_point = compute_umbilic_point(mu);
                RealVector E = point_on_side[0];
                RealVector B = point_on_side[1];
                RealVector D = point_on_side[2];

                // Add O-B as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->O(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, B, n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add W-B-G as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), B, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(B, corey->G(), n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add E-D as one Curve.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(E, D, n, temp_curve);
                    c.push_back(temp_curve);
                }  

                break;
            }
            case COREYQUAD_D_POINT:
            {
                // This stays here (and not at the ctor) because the user may have change the parameters.
                //
                std::vector<RealVector> vertex(3), point_on_side(3);
                for (int i = 0; i < 3; i++){
                    vertex_and_side(side_opposite_vertex[i], mu, vertex[i], point_on_side[i]);
                }

                RealVector umbilic_point = compute_umbilic_point(mu);
                RealVector E = point_on_side[0];
                RealVector B = point_on_side[1];
                RealVector D = point_on_side[2];

                // Add G-D as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->G(), umbilic_point, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(umbilic_point, D, n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add W-D-O as two Curves.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(corey->W(), D, n, temp_curve);
                    c.push_back(temp_curve);

                    Utilities::regularly_sampled_segment(D, corey->O(), n, temp_curve);
                    c.push_back(temp_curve);
                }  

                // Add E-B as one Curve.
                //
                {
                    Curve temp_curve;
                    Utilities::regularly_sampled_segment(E, B, n, temp_curve);
                    c.push_back(temp_curve);
                }  

                break;
            }
            default:
            {
                break;
            }
        }
    }

    return;
}

RealVector CoreyQuadExplicitHugoniotCurve::generic(void *obj, double phi){
    CoreyQuadExplicitHugoniotCurve *cqeh = (CoreyQuadExplicitHugoniotCurve*)obj;

    double muw = cqeh->muw;
    double mug = cqeh->mug;
    double muo = cqeh->muo;

    double c = cos(phi);
    double s = sin(phi);

    double c2 = c*c;
    double s2 = s*s;

    double cps = c + s;
    double cps2 = cps*cps;

    double sw = cqeh->reference_point.point(0);
    double so = cqeh->reference_point.point(1);
    double sg = 1.0 - (sw + so);

    double kw = sw*sw/muw;
    double ko = so*so/muo;
    double kg = sg*sg/mug;

    double numer1 = (sw*c/muw)*(ko + kg);
    double numer2 = kw*(so*s/muo - sg*(s + c)/mug);
    double numer3 = (so*s/muo)*(kw + kg);
    double numer4 = ko*(sw*c/muw - sg*(s + c)/mug);
    double numer  = s*(numer1 - numer2) - c*(numer3 - numer4);
           numer  = -2.*numer;

    double denom1 = (c2/muw)*(ko + kg);
    double denom2 = kw*(s2/muo + cps2/mug);
    double denom3 = (s2/muo)*(kw + kg);
    double denom4 = ko*(c2/muw + cps2/mug);
    double denom  = s*(denom1 - denom2) - c*(denom3 - denom4);

    double r;
    if (std::abs(denom) <= 1.0e-3*std::abs(numer)){
        r = 1.0e3*sign(denom)*sign(numer);
    }
    else {
        r = numer/denom;
    }

    RealVector hugxy(2);

    hugxy(0) = cqeh->reference_point.point(0) + r*c;
    hugxy(1) = cqeh->reference_point.point(1) + r*s;

    return hugxy;
}

RealVector CoreyQuadExplicitHugoniotCurve::G_bif(void *obj, double phi){
    CoreyQuadExplicitHugoniotCurve *cqeh = (CoreyQuadExplicitHugoniotCurve*)obj;

    double muw = cqeh->muw;
    double mug = cqeh->mug;
    double muo = cqeh->muo;

    double alf = cqeh->reference_point.point(0)/muw;
    double alf2 = alf*alf;

    double d = alf2*(muw + muo) + (1.- alf*(muw + muo))*(1.- alf*(muw + muo))/mug;

    double c = std::cos(phi);
    double s = std::sin(phi);

    double numer = -2.*alf2*(c + s)*(alf*(muw + muo + mug) - 1.)/mug;
    double denom = d*c*s/(muw*muo) + (c*c/muw + s*s/muo + (c + s)*(c + s)/mug )*alf2;

    double r;
    if (std::abs(denom) <= 1.0e-3*std::abs(numer)){
        r = 1.0e3*sign(denom)*sign(numer);
    }
    else {
        r = numer/denom;
    }

    RealVector hugxy(2);

    hugxy(0) = cqeh->reference_point.point(0) + r*c;
    hugxy(1) = cqeh->reference_point.point(1) + r*s;

    return hugxy;
}

RealVector CoreyQuadExplicitHugoniotCurve::O_bif(void *obj, double phi){
    CoreyQuadExplicitHugoniotCurve *cqeh = (CoreyQuadExplicitHugoniotCurve*)obj;

    double muw = cqeh->muw;
    double mug = cqeh->mug;
    double muo = cqeh->muo;

    double alf = cqeh->reference_point.point(0)/muw;
    double alf2 = alf*alf;

    double d = alf2*(muw + mug) + (1. - alf*(muw + mug))*(1. - alf*(muw + mug))/muo;

    double c = std::cos(phi);
    double s = std::sin(phi);

    double numer = -2.*alf2*(c + s)*(alf*(muw + muo + mug) -1.)/muo;
    double denom = d*c*s/(muw*mug) + (c*c/muw + s*s/mug + (c + s)*(c + s)/muo)*alf2;

    double r;
    if (std::abs(denom) <= 1.0e-3*std::abs(numer)){
        r = 1.0e3*sign(denom)*sign(numer);
    }
    else {
        r = numer/denom;
    }

    RealVector hugxy(2);

    hugxy(0) = cqeh->reference_point.point(0) + r*c;
    hugxy(1) = 1. - (cqeh->reference_point.point(0) + cqeh->reference_point.point(1)) + r*s;
    hugxy(1) = 1. - (hugxy(0) + hugxy(1));

    return hugxy;
}

RealVector CoreyQuadExplicitHugoniotCurve::W_bif(void *obj, double phi){
    CoreyQuadExplicitHugoniotCurve *cqeh = (CoreyQuadExplicitHugoniotCurve*)obj;

    double muw = cqeh->muw;
    double mug = cqeh->mug;
    double muo = cqeh->muo;

    double alf = cqeh->reference_point.point(1)/muo;
    double alf2 = alf*alf;

    double d = alf2*(mug + muo) + (1. - alf*(mug + muo) )*(1. - alf*(mug + muo) )/muw;

    double c = std::cos(phi);
    double s = std::sin(phi);

    double numer = -2.*alf2*(c + s)*(alf*(muw + muo + mug) - 1.)/muw;
    double denom = d*c*s/(mug*muo) + (c*c/mug + s*s/muo + (c + s)*(c + s)/muw)*alf2;

    double r;
    if (std::abs(denom) <= 1.0e-3*std::abs(numer)){
        r = 1.0e3*sign(denom)*sign(numer);
    }
    else {
        r = numer/denom;
    }

    RealVector hugxy(2);

    hugxy(0) = 1. - (cqeh->reference_point.point(0) + cqeh->reference_point.point(1)) + r*c;
    hugxy(1) = cqeh->reference_point.point(1) + r*s;
    hugxy(0) = 1. - (hugxy(0) + hugxy(1));

    return hugxy;
}

RealVector CoreyQuadExplicitHugoniotCurve::GO_side(void *obj, double phi){
    CoreyQuadExplicitHugoniotCurve *cqeh = (CoreyQuadExplicitHugoniotCurve*)obj;

    double muw = cqeh->muw;
    double mug = cqeh->mug;
    double muo = cqeh->muo;

//    double sw_ref = 0.0;
    double so_ref = cqeh->reference_point.point(1); // r in Cido's text.

    double so_ref2 = so_ref*so_ref;

    double D_so_ref = so_ref2/muo + (1.0 - so_ref)*(1.0 - so_ref)/mug;
    double a = (mug + muw)*so_ref2;
    double b = muw*(mug + muo)*so_ref2/muo - mug*muw*D_so_ref;

    double f = muw*so_ref2;
    double e = 2.0*f;

    double c = e + mug*muo*D_so_ref;
    double d = e + mug*muo*so_ref*D_so_ref;

    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    double A = a*cos_phi*cos_phi + b*sin_phi*sin_phi + c*cos_phi*sin_phi;
    double B = (2.0*b*so_ref - e)*sin_phi + (c*so_ref - d)*cos_phi;

    double rho = -B/A;

    RealVector hugxy(2);

    hugxy(0) = rho*cos_phi;
    hugxy(1) = so_ref + rho*sin_phi;

    return hugxy;
}

RealVector CoreyQuadExplicitHugoniotCurve::GW_side(void *obj, double phi){
    CoreyQuadExplicitHugoniotCurve *cqeh = (CoreyQuadExplicitHugoniotCurve*)obj;

    double muw = cqeh->muw;
    double mug = cqeh->mug;
    double muo = cqeh->muo;

    double sw_ref = cqeh->reference_point.point(0); // r in Cido's text.
    double so_ref = cqeh->reference_point.point(1);
    double r = 1.0 - sw_ref - so_ref;

    double r2 = r*r;

    double D_r = r2/mug + (1.0 - r)*(1.0 - r)/muw;
    double a = (muo + muw)*r2;
    double b = muw*(mug + muo)*r2/mug - muo*muw*D_r;

    double f = muo*r2;
    double e = 2.0*f;

    double c = e + mug*muw*D_r;
    double d = e + mug*muw*r*D_r;

    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    double A = a*cos_phi*cos_phi + b*sin_phi*sin_phi + c*cos_phi*sin_phi;
    double B = (2.0*b*r - e)*sin_phi + (c*r - d)*cos_phi;

    double rho = -B/A;

    RealVector hugxy(2);

    hugxy(1) = rho*cos_phi;
    hugxy(0) = 1.0 - (r + rho*(cos_phi + sin_phi));

    return hugxy;
}

RealVector CoreyQuadExplicitHugoniotCurve::WO_side(void *obj, double phi){
    CoreyQuadExplicitHugoniotCurve *cqeh = (CoreyQuadExplicitHugoniotCurve*)obj;

    double muw = cqeh->muw;
    double mug = cqeh->mug;
    double muo = cqeh->muo;

    double r = cqeh->reference_point.point(0); // r in Cido's text.

    double r2 = r*r;

    double D_r = r2/muw + (1.0 - r)*(1.0 - r)/muo;
    double a = (muo + mug)*r2;
    double b = mug*(muw + muo)*r2/muw - muo*mug*D_r;

    double f = mug*r2;
    double e = 2.0*f;

    double c = e + muo*muw*D_r;
    double d = e + muo*muw*r*D_r;

    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    double A = a*cos_phi*cos_phi + b*sin_phi*sin_phi + c*cos_phi*sin_phi;
    double B = (2.0*b*r - e)*sin_phi + (c*r - d)*cos_phi;

    double rho = -B/A;

    RealVector hugxy(2);

    hugxy(0) = r + rho*sin_phi;
    hugxy(1) = 1.0 - hugxy(0) - rho*cos_phi;

    return hugxy;
}

bool CoreyQuadExplicitHugoniotCurve::f_asymptote(void *obj, const RealVector &p, const RealVector &q){
    CoreyQuadExplicitHugoniotCurve *cqehc = (CoreyQuadExplicitHugoniotCurve*)obj;

    RealVector p_minus_ref = p - cqehc->reference_point.point;
    RealVector q_minus_ref = q - cqehc->reference_point.point;

    double prod = norm(p_minus_ref)*norm(q_minus_ref);

    return std::abs(prod + p_minus_ref*q_minus_ref) < std::abs(prod)*1e-1;
}

void CoreyQuadExplicitHugoniotCurve::vertex_and_side(int side_opposite_vertex, const RealVector &mu, RealVector &vertex, RealVector &point_on_side){
    vertex.resize(2);
    point_on_side.resize(2);

    double muw = mu(0);
    double muo = mu(1);
    double mug = mu(2);

    if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SW_ZERO){
        vertex(0) = 1.0;
        vertex(1) = 0.0;

        point_on_side(0) = 0.0;
        point_on_side(1) = muo/(mug + muo); // Equation of line from            sg/mug = so/muo
                                            //                       (1 - sw - so)/mug = so/muo, or, (1 - sw - so)*muo - so*mug = 0.
    }
    else if (side_opposite_vertex == THREE_PHASE_BOUNDARY_SO_ZERO){
        vertex(0) = 0.0;
        vertex(1) = 1.0;

        point_on_side(0) = muw/(muw + mug); // Equation of line from            sw/muw = sg/mug,
                                            //                                  sw/muw = (1 - sw - so)/mug, or, sw*mug - (1 - sw - so)*muw = 0.
        point_on_side(1) = 0.0;
    }
    else {
        vertex(0) = 0.0;
        vertex(1) = 0.0;

        point_on_side(0) = muw/(muw + muo); // Equation of line from            sw/muw = so/muo, or, sw*muo - so*muw = 0.
        point_on_side(1) = muo/(muw + muo);
    }  

    return;
}

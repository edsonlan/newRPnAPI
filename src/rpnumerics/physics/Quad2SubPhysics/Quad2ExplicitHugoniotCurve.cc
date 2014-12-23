#include "Quad2ExplicitHugoniotCurve.h"
#include "Quad2SubPhysics.h"

Quad2ExplicitHugoniotCurve::Quad2ExplicitHugoniotCurve(Quad2SubPhysics *q) : HugoniotCurve(q->flux(), q->accumulation(), q->boundary()), quad2(q) {
    method_ = EXPLICIT_HUGONIOT;

    info_ = std::string("Quad2 Explicit");
}

Quad2ExplicitHugoniotCurve::~Quad2ExplicitHugoniotCurve(){
}

void Quad2ExplicitHugoniotCurve::curve(const ReferencePoint &ref, int type, std::vector<Curve> &c){
    c.clear();

    reference_point = ref;

    int n = 100;

    double phi_begin = 0.0; // + M_PI*.1; // Ugliest possible hack to avoid problems between ParametricPlot and RectBoundary.
    double phi_end   = M_PI; // + M_PI*.1;

    // Some parameters.
    //
    alpha0 = .5*(quad2->a2() - quad2->b1());
    beta0  = .5*(quad2->b2() - quad2->c1());
    gamma0 = .5*(quad2->d2() - quad2->e1());
    alpha1 = .5*(quad2->a2() + quad2->b1());
    beta1  = .5*(quad2->b2() + quad2->c1());
    gamma1 = .5*(quad2->d2() + quad2->e1());
    alpha2 = .5*(quad2->b2() - quad2->a1());
    beta2  = .5*(quad2->c2() - quad2->b1());
    gamma2 = .5*(quad2->e2() - quad2->d1());

    std::vector<Curve> curve;
    ParametricPlot::plot(&generic, 0 /*&f_asymptote*/, (void*)this, phi_begin, phi_end, n, boundary, curve);

    for (int i = 0; i < curve.size(); i++){
        // When the asymptotes are available, the test below will be, hopefully, not needed.
        //
        Curve temp;
        int pos = 0;
        bool asymptote_found = false;

        // Skip the point if an asymptote is found.
        //
        while (pos < curve[i].curve.size() - 1){
            if (norm(curve[i].curve[pos] - curve[i].curve[pos + 1]) < 1e-1){
                temp.curve.push_back(curve[i].curve[pos]);
                asymptote_found = false;
            }
            else {
                if (temp.curve.size() > 1) c.push_back(temp);
                asymptote_found = true;
                temp.curve.clear();
            }

            pos++;
        }

        // Add the last point of the curve, which was not tested due to the 
        //
        //     curve[i].curve.size() - 1
        // 
        // above.
        //
        if (!asymptote_found){
            temp.curve.push_back(curve[i].curve.back());
            c.push_back(temp);
        }
    }

    return;
}

RealVector Quad2ExplicitHugoniotCurve::generic(void *obj, double phi){
    Quad2ExplicitHugoniotCurve *q2eh = (Quad2ExplicitHugoniotCurve*)obj;

    double cosphi = std::cos(phi);
    double sinphi = std::sin(phi);

    double c2phi = cosphi*cosphi - sinphi*sinphi;
    double s2phi = 2.*cosphi*sinphi;

    double alpha = q2eh->alpha1*c2phi + q2eh->alpha2*s2phi + q2eh->alpha0;
    double beta  = q2eh->beta1*c2phi  + q2eh->beta2*s2phi  + q2eh->beta0;
    double gamma = q2eh->gamma1*c2phi + q2eh->gamma2*s2phi + q2eh->gamma0;

    double ul = q2eh->reference_point.point(0);
    double vl = q2eh->reference_point.point(1);

    double denom = alpha*cosphi + beta*sinphi;
    double numer =  -2.0*(alpha*ul + beta*vl + gamma);

    double r;
    if (std::abs(denom) <= 1.0e-3*std::abs(numer)){
        r = 1.0e3*sign(denom)*sign(numer);
    }
    else {
        r = numer/denom;
    }

    RealVector hugxy(2);

    hugxy(0) = q2eh->reference_point.point(0) + r*cosphi;
    hugxy(1) = q2eh->reference_point.point(1) + r*sinphi;

    return hugxy;
}

//bool Quad2ExplicitHugoniotCurve::f_asymptote(void *obj, const RealVector &p, const RealVector &q){
//    Quad2ExplicitHugoniotCurve *cqehc = (Quad2ExplicitHugoniotCurve*)obj;

//    RealVector p_minus_ref = p - cqehc->reference_point.point;
//    RealVector q_minus_ref = q - cqehc->reference_point.point;

//    double prod = norm(p_minus_ref)*norm(q_minus_ref);

//    return std::abs(prod + p_minus_ref*q_minus_ref) < std::abs(prod)*1e-1;
//}


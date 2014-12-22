#include "ClassifyShockCurve.h"

ClassifyShockCurve::ClassifyShockCurve(const HugoniotContinuation *h) : hc(h) {
    sp = std::string("+"); // sigma > lambda
    sm = std::string("-"); // sigma < lambda
    sd = std::string("."); // lambda is complex
}

ClassifyShockCurve::~ClassifyShockCurve(){
}

// This method is the engine of classify_point(). Classify based on a given set of eigenpairs.
//
int ClassifyShockCurve::half_classify(double sigma, const std::vector<eigenpair> &e, std::string &s){
    for (int i = 0; i < e.size(); i++){
        // If this eigenvalue is complex, the next one will also be complex. 
        // One is the conjugate of the other.
        if (fabs(e[i].i) > 1e-10){
            s += sd; // These two are is complex
            if (sigma > e[i].r) s += sp;
            else                s += sm;

            // Skip the next one.
            i++;
        }
        else {
            if (sigma > e[i].r) s += sp;
            else                s += sm;
        }
    }

    // TODO: Add a meaningful return code.
    return 1;
}

int ClassifyShockCurve::classify_point(const RealVector &p, const ReferencePoint &ref, std::string &s){
    int n = p.size();

    WaveState u(p);

    JetMatrix jet_F(n), jet_G(n);

    hc->flux()->jet(u, jet_F, 1);
    hc->accumulation()->jet(u, jet_G, 1);

    // Obtain sigma
    double sigma = hc->sigma(jet_F.function(), jet_G.function());

    // Obtain the lambdas at the current point
    std::vector<eigenpair> e;
    Eigen::eig(n, jet_F.Jacobian().data(), jet_G.Jacobian().data(), e);

    // Classify proper
    s.clear();

    // First compare sigma with lambda_ref...
    int info_ref = half_classify(sigma, ref.e, s);

    // ...then with the lambdas at the point.
    int info_current = half_classify(sigma, e, s);

    return 0;
}


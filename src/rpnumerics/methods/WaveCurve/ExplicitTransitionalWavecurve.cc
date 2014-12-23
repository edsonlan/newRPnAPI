#include "ExplicitTransitionalWavecurve.h"

ExplicitTransitionalWavecurve::ExplicitTransitionalWavecurve(){
    corey_ = new CoreyQuadSubPhysics;
}

ExplicitTransitionalWavecurve::~ExplicitTransitionalWavecurve(){
    delete corey_;
}

void ExplicitTransitionalWavecurve::transitional_double_contact_jet(double L, int degree, JetMatrix &Ljet){
    Ljet.resize(1);

    if (degree >= 0){
        double L2 = L*L;

        double r2 = r*r;
        double r3 = r2*r;

        double func;

//        double func = (1.0 - L)*(r2*L2 + r*pow(2.0*(1.0 + r)*L2 - 3.0*r*L + r, 2.0)) -
//                   (2.0*(1.0 + r)*L2 - 2.0*r*L + r)*pow(L2 + r*pow(1.0 - L, 2.0), 2.0);

//        std::cout << "Polynomial. r = " << r << ", L = " << L << ", pol. = " << func << std::endl;

        double a = -(6.0*r2 + 2.0*r3 + 2.0 + 6.0*r);
        double b = 6.0*r + 12.0*r2 + 6.0*r3;
        double c = -(r + 6.0*r2 + 5.0*r3);
        double d = -r2 - r3;
        double e = r2 + 3.0*r3;
        double f = -r3;

        func = f + L*(e + L*(d + L*(c + L*(b + L*a))));

//        std::cout << "Second version, pol. = " << func << std::endl;

        Ljet.set(0, func);

        if (degree >= 1){
            double df = e + L*(2.0*d + L*(3.0*c + L*(4.0*b + L*5.0*a)));
            Ljet.set(0, 0, df);
        }
    }

    return;
}

int ExplicitTransitionalWavecurve::L_function(void *obj, double L, double &s_C2){
    ExplicitTransitionalWavecurve *etwc = (ExplicitTransitionalWavecurve*)obj;

    JetMatrix Ljet;
    etwc->transitional_double_contact_jet(L, 0, Ljet);

    s_C2 = Ljet.get(0);

    return BISECTION_FUNCTION_OK;
}

void ExplicitTransitionalWavecurve::U_side_polynomial_jet(double L, int degree, JetMatrix &Ljet){
    Ljet.resize(1);

    if (degree >= 0){
        double L2 = L*L;
        double L3 = L2*L;

        double p = (1.0 + r)*(1.0 - 2.0*sM)*L3 + (1.0 + r)*sM*(1.0 + 2.0*sM)*L2 - (r + 2.0*(1.0 + r)*sM*sM)*L + r*sM;

        Ljet.set(0, p);

        if (degree >= 1){
            double dp = 0.0; // TODO: Compute this.
            Ljet.set(0, 0, dp);
        }
    }

    return;
}

int ExplicitTransitionalWavecurve::U_side_function(void *obj, double L, double &sl){
    ExplicitTransitionalWavecurve *etwc = (ExplicitTransitionalWavecurve*)obj;

    JetMatrix Ljet;
    etwc->U_side_polynomial_jet(L, 0, Ljet);

    sl = Ljet.get(0);

    return BISECTION_FUNCTION_OK;
}

RealVector ExplicitTransitionalWavecurve::EW_line(double sM){
    RealVector p(2);

    double muo = corey_->muo()->value();
    double mug = corey_->mug()->value(); 

    p(0) = 1.0 - sM;
    p(1) = muo*sM/(muo + mug);

    return p;
}

double ExplicitTransitionalWavecurve::EW_s(const RealVector &M){
    return 1.0 - M(0);
}

RealVector ExplicitTransitionalWavecurve::BO_line(double sM){
    RealVector p(2);

    double muw = corey_->muw()->value();
    double mug = corey_->mug()->value();

    p(0) = muw*sM/(muw + mug);
    p(1) = 1.0 - sM;

    return p;
}

double ExplicitTransitionalWavecurve::BO_s(const RealVector &M){
    return 1.0 - M(1);
}

RealVector ExplicitTransitionalWavecurve::DG_line(double sM){
    RealVector p(2);

    double muw = corey_->muw()->value();
    double muo = corey_->muo()->value();

    p(0) = muw*sM/(muw + muo);
    p(1) = muo*sM/(muw + muo);

    return p;
}

double ExplicitTransitionalWavecurve::DG_s(const RealVector &M){
    return M(0) + M(1);
}

// This version is deprecated.
//
//int ExplicitTransitionalWavecurve::subdivide_curve(int type, std::vector<double> &s_param, std::vector<RealVector> &extpts){
//    s_param.clear();
//    extpts.clear();

//    double muw = corey_->muw()->value();
//    double muo = corey_->muo()->value();
//    double mug = corey_->mug()->value();

//    double mu_total = muw + muo + mug;

//    // Vertex (all cases are treated equally).
//    //
//    double s_vertex = 0.0;

//    // E2 and E1.
//    //
//    double s_E2, s_E1;

//    if (type == EW){
//        s_E1 = (muw + mug)/(mu_total + muw);
//        s_E2 = 1.0 - sqrt(muw/mu_total);

//        // r (class member).
//        //
//        r = (mug + muo)/muw;
//    }
//    else if (type == DG){
//        s_E1 = (muo + mug)/(mu_total + mug);
//        s_E2 = 1.0 - sqrt(mug/mu_total);

//        // r (class member).
//        //
//        r = (muw + muo)/mug;
//    }
//    else if (type == BO){
//        s_E1 = (muw + muo)/(mu_total + muo);
//        s_E2 = 1.0 - sqrt(muo/mu_total);

//        // r (class member).
//        //
//        r = (muw + mug)/muo;
//    }
//    
//    // C2 and C1.
//    // Find the only root that is supposed to exist in the interval [0, 1].
//    //
//    double s_min = 0.0;
//    double s_max = 1.0;

//    int n = 100;
//    double delta = (s_max - s_min)/(double)(n - 1);

//    double L_prev;
//    double s_prev = s_min;

//    int info_prev = L_function((void*)this, s_prev, L_prev);
//    if (info_prev == BISECTION_FUNCTION_ERROR) return CURVE_ERROR;

//    int pos = 1;
//    bool found = false;

//    while (pos < n && !found){
//        double L;
//        double s = s_min + delta*(double)pos;

//        int info_L = L_function((void*)this, s, L);
//        if (info_L == BISECTION_FUNCTION_ERROR) return CURVE_ERROR;

//        if (L*L_prev < 0.0){
//            s_min = s_prev;
//            s_max = s;

//            found = true;
//        }
//        else {
//            L_prev = L;
//            s_prev = s;

//            pos++;
//        }
//    }

//    if (!found) return CURVE_ERROR;

//    // Cido is trying to find the roots of the polynomial. If he succeeds then the step below
//    // can be safely replaced by the appropriate root. If he doesn't, the Newton method
//    // can replace the bisection below, because now the derivative of the polynomial
//    // is also known.
//    //
//    double s_C2;
//    int info_L = Bisection::simple_bisection(s_min, s_max, (void*)this, &L_function, s_C2);

//    double s_C1 = r*s_C2/(2.0*(1.0 + r)*s_C2*s_C2 - 2.0*r*s_C2 + r);

//    double s_umbilic;
//    if      (type == EW) s_umbilic = (mug + muo)/mu_total;
//    else if (type == DG) s_umbilic = (muo + muw)/mu_total;
//    else if (type == BO) s_umbilic = (mug + muw)/mu_total;

//    // Point on side opposing vertex (all cases are treated equally).
//    //
//    double s_side = 1.0;

//    s_param.push_back(s_vertex);
//    s_param.push_back(s_E2);
//    s_param.push_back(s_E1);
//    s_param.push_back(s_C2);
//    s_param.push_back(s_umbilic);
//    s_param.push_back(s_C1);
//    s_param.push_back(s_side);

//    if (type == EW){
//        for (int i = 0; i < s_param.size(); i++) extpts.push_back(EW_line(s_param[i]));
//    }
//    else if (type == BO){
//        for (int i = 0; i < s_param.size(); i++) extpts.push_back(BO_line(s_param[i]));
//    }
//    if (type == DG){
//        for (int i = 0; i < s_param.size(); i++) extpts.push_back(DG_line(s_param[i]));
//    }

//    return CURVE_OK;
//}

// The names below match those in Cido & Freddie's paper.
//
// Vector extpts is of the form:
// 
//     extpts[0] = Vertex (W, O or G).
//     extpts[1] = E2.
//     extpts[2] = E1.
//     extpts[3] = C1.
//     extpts[4] = Umbilic point.
//     extpts[5] = C2.
//     extpts[6] = Point on the side opposing the vertex (E, B or D).
//
// This method should be cleaned up, by abstracting the mu's.
//
int ExplicitTransitionalWavecurve::subdivide_curve(int type, std::vector<double> &s_param, std::vector<RealVector> &extpts){
    s_param.clear();
    extpts.clear();

    double muw = corey_->muw()->value();
    double muo = corey_->muo()->value();
    double mug = corey_->mug()->value();

    double mu_total = muw + muo + mug;

    // Vertex (all cases are treated equally).
    //
    double s_vertex = 0.0;

    // E2 and E1.
    //
    double s_E2, s_E1;

    if (type == EW){
        s_E1 = (muw + mug)/(mu_total + muw);
        s_E2 = 1.0 - sqrt(muw/mu_total);

        // r (class member).
        //
        r = (mug + muo)/muw;
    }
    else if (type == DG){
        s_E1 = (muo + mug)/(mu_total + mug);
        s_E2 = 1.0 - sqrt(mug/mu_total);

        // r (class member).
        //
        r = (muw + muo)/mug;
    }
    else if (type == BO){
        s_E1 = (muw + muo)/(mu_total + muo);
        s_E2 = 1.0 - sqrt(muo/mu_total);

        // r (class member).
        //
        r = (muw + mug)/muo;
    }
    
    // C2 and C1.
    // Find the only root that is supposed to exist in the interval [0, 1].
    //
    double root = sqrt(2.0*r/(1.0 + r));

    double s_C2 = root*.5;

    double s_C1 = -.5*root/(root - 2.0);

    double s_umbilic;
    if      (type == EW) s_umbilic = (mug + muo)/mu_total;
    else if (type == DG) s_umbilic = (muo + muw)/mu_total;
    else if (type == BO) s_umbilic = (mug + muw)/mu_total;

    // Point on side opposing vertex (all cases are treated equally).
    //
    double s_side = 1.0;

    s_param.push_back(s_vertex);
    s_param.push_back(s_E2);
    s_param.push_back(s_E1);
    s_param.push_back(s_C2);
    s_param.push_back(s_umbilic);
    s_param.push_back(s_C1);
    s_param.push_back(s_side);

    if (type == EW){
        for (int i = 0; i < s_param.size(); i++) extpts.push_back(EW_line(s_param[i]));
    }
    else if (type == BO){
        for (int i = 0; i < s_param.size(); i++) extpts.push_back(BO_line(s_param[i]));
    }
    if (type == DG){
        for (int i = 0; i < s_param.size(); i++) extpts.push_back(DG_line(s_param[i]));
    }

    return CURVE_OK;
}

int ExplicitTransitionalWavecurve::find_subdivision(int type, const RealVector &M, const std::vector<double> &s_param, RealVector &point_sl, RealVector &point_sr){
    double muw = corey_->muw()->value();
    double muo = corey_->muo()->value();
    double mug = corey_->mug()->value();

    double mu_total = muw + muo + mug;

    if (type == EW){
        // Class member.
        //
        sM = EW_s(M);

        // Class member.
        //
        r = (mug + muo)/muw;
    }
    else if (type == BO){
        sM = BO_s(M);

        // r (class member).
        //
        r = (muw + mug)/muo;
    }
    else if (type == DG){
        sM = DG_s(M);

        // r (class member).
        //
        r = (muw + muo)/mug;
    }

    int pos = 0;
    bool found = false;
    while (pos < s_param.size() - 1 && !found){
        if (sM >= s_param[pos] && sM < s_param[pos + 1]) found = true;
        else pos++;
    }

    if (!found) return SUBDIVISION_ERROR;

    double sl, sr;
    bool interval_found = false;
    double s_umbilic = s_param[4];

    if (pos == 0){
        // M lies between V and E2.
        // Nothing to be done here: no interval can be found.
    }
    else if (pos == 1){
        std::cout << "M lies between E2 and E1." << std::endl;

        // M lies between E2 and E1.
        //
        double sigma = 2.0*sM*(1.0 - sM)/(r*pow(sM*sM/r + pow(1.0 - sM, 2.0), 2.0));
        double aux = sigma*(pow(sM, 2.0) + r*pow(1.0 - sM, 2.0));

        double A = -aux*(1.0 + r);
        double B = aux*2.0*r + r*(1.0 - 2.0*sM);
        double C = -aux*r + r*sM;

        double x1(3.14), x2(3.14);
        int info = Utilities::Bhaskara(B/A, C/A, x1, x2);
        std::cout << "    E1--E2. A = "<< A << ", B = " << B << ", C = " << C << std::endl;
        std::cout << "    E1--E2. x1 = "<< x1 << ", x2 = " << x2 << ", info = " << info << std::endl;

        if (info != BHASKARA_TWO_DIFFERENT_ROOTS) return SUBDIVISION_ERROR;

        std::cout << "s_umbilic = " << s_umbilic << std::endl;

        double s_umbilic_line = s_umbilic;
//        s_umbilic_line = 0.0; // Hack.

        if ((x1 >= s_umbilic_line && x1 <= 1.0) &&
            (x2 >= s_umbilic_line && x2 <= 1.0)) return SUBDIVISION_ERROR;

        if      (x1 >= s_umbilic_line && x1 <= 1.0) sl = x1;
        else if (x2 >= s_umbilic_line && x2 <= 1.0) sl = x2;
        else return SUBDIVISION_ERROR;

        sr = 1.0;

        interval_found = true;

//        // Fast Lambda @ M:
//        //
//        const FluxFunction *f = corey_->flux();
//        JetMatrix fMjet(2);
//        f->jet(M, fMjet, 1);

//        std::vector<eigenpair> e;
//        Eigen::eig(2, fMjet.Jacobian().data(), e);

//        double fast_lambda_M = e[1].r;

//        // Find the point where sigma == fast_lambda_M.
//        //
//        double s_min = 0.0, s_max = 1.0;
//        int n = 100;
//        int j = 0;
//        double delta = (s_max - s_min)/(int)(n - 1);
//        bool found = false;
//        double lambda_minus_sigma_prev;

//        while (j < n && !found){
//            double s_sigma = s_min + delta*(double)j;

//            RealVector candidate;
//            if      (type == EW) candidate = EW_line(s_sigma);
//            else if (type == DG) candidate = DG_line(s_sigma);
//            else if (type == BO) candidate = BO_line(s_sigma);

//            JetMatrix Cjet(2);
//            f->jet(candidate, Cjet, 0);

//            double den = 0.0;
//            double num = 0.0;

//            for (int i = 0; i < 2; i++){
//                den += (candidate(i) - M(i))*(candidate(i) - M(i));
//                num += (candidate(i) - M(i))*(Cjet.function()(i) - fMjet.function()(i));
//            }

//            double sigma = num/den;
//            double lambda_minus_sigma = fast_lambda_M - sigma;

//            if (j > 0){
//                if (lambda_minus_sigma_prev*lambda_minus_sigma < 0.0) {
//                    std::cout << "Found! s_sigma = " << s_sigma << std::endl;
//                }
//            }

//            j++;

//            lambda_minus_sigma_prev = lambda_minus_sigma;
//        }

    }
    else if (pos == 2){
        // M lies between E1 and C2.
        //
        double sigma = 2.0*sM*(1.0 - sM)/(r*pow(sM*sM/r + pow(1.0 - sM, 2.0), 2.0));
        double aux = sigma*(pow(sM, 2.0) + r*pow(1.0 - sM, 2.0));

        double A = -aux*(1.0 + r);
        double B = aux*2.0*r + r*(1.0 - 2.0*sM);
        double C = -aux*r + r*sM;

        double x1, x2;
        int info = Utilities::Bhaskara(B/A, C/A, x1, x2);
        std::cout << "E1--C2. sl, x1 = " << x1 << ", x2 = " << x2 << std::endl;

        if (info != BHASKARA_TWO_DIFFERENT_ROOTS) return SUBDIVISION_ERROR;

        double s_C1 = s_param[5];

        if ((x1 >= s_C1 && x1 <= 1.0) &&
            (x2 >= s_C1 && x2 <= 1.0)) return SUBDIVISION_ERROR;

        if      (x1 >= s_C1 && x1 <= 1.0) sl = x1;
        else if (x2 >= s_C1 && x2 <= 1.0) sl = x2;
        else return SUBDIVISION_ERROR;

        std::cout << "E1--C2. sl = " << sl << std::endl;

        // sr.
        //
        sigma = 2.0*sM/(sM*sM + r*pow(1.0 - sM, 2.0));
        aux = sigma*(pow(sM, 2.0) + r*pow(1.0 - sM, 2.0));
        A = -aux*(1.0 + r);
        B = aux*2.0*r + r*(1.0 - 2.0*sM);
        C = -aux*r + r*sM;

        info = Utilities::Bhaskara(B/A, C/A, x1, x2);

        if (info != BHASKARA_TWO_DIFFERENT_ROOTS) return SUBDIVISION_ERROR;

        if ((x1 >= sl && x1 <= 1.0) &&
            (x2 >= sl && x2 <= 1.0)) return SUBDIVISION_ERROR;

        if      (x1 >= sl && x1 <= 1.0) sr = x1;
        else if (x2 >= sl && x2 <= 1.0) sr = x2;
        else return SUBDIVISION_ERROR;

        interval_found = true;
    }
    else if (pos == 3){
        // M lies between C2 and U.
        //
        double A = -2.0*(sM*sM + r*pow(1.0 - sM, 2.0)) + r - 2.0*r*sM;
        double B = r*sM;

        sl = -B/A;

        std::cout << "C2--U. sl = " << sl << std::endl;

        if      (sl < 0.0) sl = 0.0;
        else if (sl > 1.0) sl = 1.0;

        double sigma = 2.0*sM/(sM*sM + r*pow(1.0 - sM, 2.0));
        double aux = sigma*(pow(sM, 2.0) + r*pow(1.0 - sM, 2.0));
        A = -aux*(1.0 + r);
        B = aux*2.0*r + r*(1.0 - 2.0*sM);
        double C = -aux*r + r*sM;
 
        double x1, x2;
        int info = Utilities::Bhaskara(B/A, C/A, x1, x2);
        if (info != BHASKARA_TWO_DIFFERENT_ROOTS) return SUBDIVISION_ERROR;

        if ((x1 >= sl && x1 <= 1.0) &&
            (x2 >= sl && x2 <= 1.0)) return SUBDIVISION_ERROR;

        if      (x1 >= sl && x1 <= 1.0) sr = x1;
        else if (x2 >= sl && x2 <= 1.0) sr = x2;
        else return SUBDIVISION_ERROR;

        interval_found = true;
    }
    else if (pos == 4 || pos == 5){
        // M lies between U and C1 or between C1 and E, B, D.
        //
        double A = -2.0*(sM*sM + r*pow(1.0 - sM, 2.0)) + r - 2.0*r*sM; 
        double B = r*sM;

        sr = -B/A; 

        std::cout << "U--E: sr = " << sr << std::endl;

        double s_min = s_param[1];
        double s_max = sr;

        // Deprecated (replaced by an exact solution, viz., raiz8.
        // CIDO, please check this!
        //
        // int info_sl = Bisection::simple_bisection(s_min, s_max, (void*)this, &U_side_function, sl);
        // if (info_sl != BISECTION_CONVERGENCE_OK) return CURVE_ERROR;
        double inv_1pr = 1.0/(1.0 + r);
//        sl = (sM - sqrt(sM*sM + r*(1.0 - 2.0*sM)*inv_1pr) )/(2.0*sM - 1.0);
        sl = r*inv_1pr / (sM + sqrt(sM*sM + r*(1.0 - 2.0*sM)*inv_1pr) ); // Suggested by Cido.
                
        std::cout << "U--E: s_min = " << s_min << ", sl = " << sl << std::endl;

        interval_found = true;
    }

    if (type == EW){
        point_sl = EW_line(sl);
        point_sr = EW_line(sr);
    }
    else if (type == BO){
        point_sl = BO_line(sl);
        point_sr = BO_line(sr);
    }
    else if (type == DG){
        point_sl = DG_line(sl);
        point_sr = DG_line(sr);
    }

    if (interval_found) return SUBDIVISION_OK;
    else                return SUBDIVISION_ERROR;
}


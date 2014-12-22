#include "FoamViscosity.h"

FoamViscosity::FoamViscosity(Parameter *mug0, 
                             Parameter *epdry,
                             Parameter *fdry,
                             Parameter *foil,
                             Parameter *fmdry,
                             Parameter *fmmob,
                             Parameter *fmoil,
                             ThreePhaseFlowSubPhysics *t):
                             ThreePhaseFlowViscosity(t),
                             mug0_parameter(mug0),
                             epdry_parameter(epdry),
                             fdry_parameter(fdry),
                             foil_parameter(foil),
                             fmdry_parameter(fmdry),
                             fmmob_parameter(fmmob),
                             fmoil_parameter(fmoil){
}

FoamViscosity::~FoamViscosity(){
}

// This polynomial which provides a Heaviside-like with
// two continuous derivatives in the [-1, 1] interval.
// The first and second derivatives are zero at the interval's
// endpoints.
//
void FoamViscosity::Fdry(double x, int degree, JetMatrix &fdry_jet){
    fdry_jet.resize(1);

    if (degree >= 0){
        double x2 = x*x;
        double x4 = x2*x2;

        fdry_jet.set(0, .375*x*(5.0 - 10.0*x2/3.0 + x4));

        if (degree >= 1){
            fdry_jet.set(0, 0, 1.875 - 3.75*x2 + 1.875*x4);

            if (degree >= 2){
                fdry_jet.set(0, 0, 0, 7.5*(x2*x - x));
            }
        }
    }

    return;
}

void FoamViscosity::Fo(double so, int degree, JetMatrix &fo_jet){
    double epoil = 2.0;
    double fmoil = 1.0;
    double floil = 0.0;

    fo_jet.resize(1);

    if (so >= fmoil){
        fo_jet.set(0, 0.0);
        fo_jet.set(0, 0, 0.0);
        fo_jet.set(0, 0, 0, 0.0);
    }
    else if (so <= floil){
        fo_jet.set(0, 1.0);
        fo_jet.set(0, 0, 0.0);
        fo_jet.set(0, 0, 0, 0.0);
    }
    else {
        double inv = 1.0/(fmoil - floil);
        double y = (fmoil - so)*inv;
        double fo2 = pow(y, epoil - 2.0);
        double fo1 = fo2*y;
        double fo  = fo2*fo2;

        if (degree >= 0){
            fo_jet.set(0, fo);

            if (degree >= 1){
                fo_jet.set(0, 0, -epoil*fo1*inv);

                if (degree >= 2){
                    fo_jet.set(0, 0, 0, epoil*(epoil - 1.0)*fo2*inv*inv);
                }
            }
        }
    }

    return;
}

int FoamViscosity::gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet){
    double epdry = epdry_parameter->value();
    double fdry  = fdry_parameter->value();
    double foil  = foil_parameter->value();
    double fmdry = fmdry_parameter->value();
    double fmmob = fmmob_parameter->value();
    double fmoil = fmoil_parameter->value();

    mug_jet.resize(2, 1);

    JetMatrix fdry_jet;
    Fdry(w(0) - fmdry, degree, fdry_jet);

    JetMatrix fo_jet;
    Fo(w(1), degree, fo_jet);

    if (degree >= 0){
        double mug0 = mug0_parameter->value();
        double fdry = fdry_jet.get(0);
        double fo   = fo_jet.get(0);

        double mug  = mug0*(1.0 + fmmob*fdry*fo); // Fs is missing.

        mug_jet.set(0, mug);

        if (degree >= 1){
            double dfdry_dsw = fdry_jet.get(0, 0);
            double dfo_dso   = fo_jet.get(0, 0);

            mug_jet.set(0, 0, mug0*fmmob*dfdry_dsw*fo);
            mug_jet.set(0, 1, mug0*fmmob*fdry*dfo_dso);

            if (degree >= 2){
                double d2fdry_dsw2 = fdry_jet.get(0, 0, 0);
                double d2fo_dso2   = fo_jet.get(0, 0, 0);

                mug_jet.set(0, 0, 0, mug0*fmmob*d2fdry_dsw2*fo);
                mug_jet.set(0, 0, 1, mug0*fmmob*dfdry_dsw*dfo_dso);
                mug_jet.set(0, 1, 0, mug0*fmmob*dfdry_dsw*dfo_dso);
                mug_jet.set(0, 1, 1, mug0*fmmob*fdry*d2fo_dso2);
            }
        }
    }

    return VISCOSITY_OK;
}

//int FoamViscosity::gas_viscosity_jet(const WaveState &w, int degree, JetMatrix &mug_jet){
//    mug_jet.resize(2, 1);

//    double sg = 1.0 - w(0) - w(1) + 1e-6;

//    if (degree >= 0){
//        double mug  = 1.0/(sg*sg);

//        mug_jet.set(0, mug);

//        if (degree >= 1){
//            mug_jet.set(0, 0, 2.0/(sg*sg*sg));
//            mug_jet.set(0, 1, 2.0/(sg*sg*sg));

//            if (degree >= 2){
//                double d2 = 6.0/(sg*sg*sg*sg);
//                mug_jet.set(0, 0, 0, d2);
//                mug_jet.set(0, 0, 1, d2);
//                mug_jet.set(0, 1, 0, d2);
//                mug_jet.set(0, 1, 1, d2);
//            }
//        }
//    }

//    return VISCOSITY_OK;
//}

double FoamViscosity::gas_viscosity(const RealVector &p){
    JetMatrix mug_jet;
    gas_viscosity_jet(WaveState(p), 0, mug_jet);

    return mug_jet.get(0);
}


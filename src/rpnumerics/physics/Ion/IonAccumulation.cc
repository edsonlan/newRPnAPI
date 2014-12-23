#include "IonAccumulation.h"

IonAccumulation::IonAccumulation(const IonAdsorption *a) : AccumulationFunction(), adsorption(a) {
}

IonAccumulation::~IonAccumulation(){
}

int IonAccumulation::jet(const WaveState &w, JetMatrix &m, int degree) const {
    JetMatrix a_jet(2);
    adsorption->jet(w, a_jet, degree);

    if (degree >= 0){
        double u = w(0);
        double v = w(1);

        double a = a_jet.get(0);

        m.set(0, u);
        m.set(1, a + u*v);

        if (degree >= 1){
            double da_du = a_jet.get(0, 0);
            double da_dv = a_jet.get(0, 1);

            m.set(0, 0, 1.0);
            m.set(0, 1, 0.0);
            m.set(1, 0, da_du + v);
            m.set(1, 1, da_dv + u);

            if (degree == 2){
                double d2a_du2  = a_jet.get(0, 0, 0);
                double d2a_dudv = a_jet.get(0, 0, 1);
                double d2a_dv2  = a_jet.get(0, 1, 1);

                m.set(0, 0, 0, 0.0);
                m.set(0, 0, 1, 0.0);
                m.set(0, 1, 0, 0.0);
                m.set(0, 1, 1, 0.0);

                m.set(1, 0, 0, d2a_du2);
                m.set(1, 0, 1, d2a_dudv + 1.0);
                m.set(1, 1, 0, d2a_dudv + 1.0);
                m.set(1, 1, 1, d2a_dv2);
            }
        }
    }

    return 2;
}

IonAccumulation * IonAccumulation::clone() const{
    return new IonAccumulation(adsorption);
}


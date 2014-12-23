#include "ThreePhaseFlowAccumulation.h"

ThreePhaseFlowAccumulation::ThreePhaseFlowAccumulation(){}

ThreePhaseFlowAccumulation::~ThreePhaseFlowAccumulation(){}

int ThreePhaseFlowAccumulation::jet(const WaveState &w, JetMatrix &m, int degree) const {
    m.resize(2);

    if (degree >= 0){
        m.set(0, w(0));
        m.set(1, w(1));

        if (degree >= 1){
            m.set(0, 0, 1.0);
            m.set(0, 1, 0.0);
            m.set(1, 0, 0.0);
            m.set(1, 1, 1.0);

            if (degree == 2){
                m.set(0, 0, 0, 0.0);
                m.set(0, 0, 1, 0.0);
                m.set(0, 1, 0, 0.0);
                m.set(0, 1, 1, 0.0);

                m.set(1, 0, 0, 0.0);
                m.set(1, 0, 1, 0.0);
                m.set(1, 1, 0, 0.0);
                m.set(1, 1, 1, 0.0);
            }
        }
    }
    return 2;
}


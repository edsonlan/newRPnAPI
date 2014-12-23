#include "Quad2AccumulationFunction.h"

Quad2AccumulationFunction::Quad2AccumulationFunction(){
}

Quad2AccumulationFunction::~Quad2AccumulationFunction(void) {}

int Quad2AccumulationFunction::jet(const WaveState &w, JetMatrix &m, int degree) const {
    m.resize(2);

    if (degree >= 0) {
        m.set(0, w(0));
        m.set(1, w(1));

        if (degree >= 1) {
            m.set(0, 0, 1.0);
            m.set(0, 1, 0.0);
            m.set(1, 0, 0.0);
            m.set(1, 1, 1.0);

            if (degree == 2) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k < 2; k++) m.set(i, j, k, 0.0);
                    }
                }
            }
        }
    }

    return 2;
}


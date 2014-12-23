#include "IonPermeability.h"

IonPermeability::IonPermeability(){
}

IonPermeability::~IonPermeability(){
}

int IonPermeability::jet(const WaveState &w, JetMatrix &p, int degree) const {
    if (degree >= 0){
        double u = w(0);

        // permw
        p.set(0, u*u);

        // permg
        p.set(1, 1.0 - 2.0*u + u*u);

        if (degree >= 1){
            p.set(0, 0, 2.0*u);
            p.set(0, 1, 0.0);

            p.set(1, 0, -2.0 + 2.0*u);
            p.set(1, 1, 0.0);

            if (degree == 2){
                p.set(0, 0, 0, 2.0);
                p.set(0, 0, 1, 0.0);
                p.set(0, 1, 1, 0.0);

                p.set(1, 0, 0, 2.0);
                p.set(1, 0, 1, 0.0);
                p.set(1, 1, 1, 0.0);
            }
        }
    }

    return 2;
}


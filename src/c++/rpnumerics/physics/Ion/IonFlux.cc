#include "IonFlux.h"

IonFlux::IonFlux(const IonRatios *r, const IonPermeability *p) : FluxFunction(), ratios(r), permeability(p){
}

IonFlux::~IonFlux(){
}

int IonFlux::jet(const WaveState &w, JetMatrix &m, int degree) const {
    JetMatrix r_jet(1);
    ratios->jet(w, r_jet, degree);

    JetMatrix p_jet(2);
    permeability->jet(w, p_jet, degree);

    if (degree >= 0){
        double u = w(0);
        double v = w(1);

        double r = r_jet.get(0);

        double permw = p_jet.get(0);
        double permg = p_jet.get(1);

        // Auxiliary value
        double inv = 1.0/(permw + r*permg);

        m.set(0, permw*inv);
        m.set(1, v*permw*inv);

        if (degree >= 1){
            double dr_dv = r_jet.get(0, 0);

            double dpermw_du = p_jet.get(0, 0);
            double dpermw_dv = p_jet.get(0, 1);
            double dpermg_du = p_jet.get(1, 0);
            double dpermg_dv = p_jet.get(1, 1);

            // Auxiliary value
            double inv2 = inv*inv;

            m.set(0, 0, -r*(dpermw_du*permg + permw*dpermg_du)*inv2);
            m.set(0, 1, -(-dpermw_dv*r*permg + permw*dr_dv*permg + permw*r*dpermg_dv)*inv2);
            m.set(1, 0, -v*r*(-dpermw_du*permg + permw*dpermg_du)*inv2);
            m.set(1, 1, (permw*permw + permw*r*permg + v*dpermw_dv*r*permg - v*permw*dr_dv*permg - v*permw*r*dpermg_dv)*inv2);

            if (degree == 2){
                double d2r_dv2  = r_jet.get(0, 0, 0);

                double d2permw_du2  = p_jet.get(0, 0, 0);
                double d2permw_dudv = p_jet.get(0, 0, 1);
                double d2permw_dv2  = p_jet.get(0, 1, 1);

                double d2permg_du2  = p_jet.get(1, 0, 0);
                double d2permg_dudv = p_jet.get(1, 0, 1);
                double d2permg_dv2  = p_jet.get(1, 1, 1);

                // Auxiliary value
                double inv3 = inv2*inv;
                 
                m.set(0, 0, 0, -r*(-d2permw_du2*permg + permw*d2permg_du2)*inv2 + 
                               2.0*r*(-dpermw_du*permg + permw*dpermg_du)*inv3*(dpermw_du+r*dpermg_du)
                     );

                m.set(0, 0, 1, -dr_dv*(-dpermw_du*permg + permw*dpermg_du)*inv -
                               r*(-d2permw_dudv*permg - dpermw_du*dpermg_dv + dpermw_dv*dpermg_du + permw*d2permg_dudv)*inv2 +
                               2.0*r*(-dpermw_du*permg + permw*dpermg_du)*inv3*(dpermw_dv + dr_dv*permg + r*dpermg_dv)
                     );

                m.set(0, 1, 0, -(-d2permw_dudv*r*permg - dpermw_dv*r*dpermg_du + dpermw_du*dr_dv*permg +
                                 permw*dr_dv*dpermg_du + dpermw_du*r*dpermg_dv + permw*r*d2permg_dudv
                                )*inv2 +
                               2.0*(-dpermw_dv*r*permg + permw*dr_dv*permg + permw*r*dpermg_dv)*inv3*(dpermw_du + r*dpermg_du)
                     );

                m.set(0, 1, 1, -(-d2permw_dv2*r*permg + permw*d2r_dv2*permg + 2.0*permw*dr_dv*dpermg_dv + permw*r*d2permg_dv2)*inv2 +
                               2.0*(-dpermw_dv*r*permg + permw*dr_dv*permg + permw*r*dpermg_dv)*inv3*(dpermw_dv + dr_dv*permg + r*dpermg_dv)
                     );




                m.set(1, 0, 0, -v*r*(-d2permw_du2*permg + permw*d2permg_du2)*inv2 + 
                               2.0*v*r*(-dpermw_du*permg + permw*dpermg_du)*inv3*(dpermw_du + r*dpermg_du)
                     );

                m.set(1, 0, 1, -r*(-dpermw_du*permg + permw*dpermg_du)*inv2 -
                               v*dr_dv*(-dpermw_du*permg + permw*dpermg_du)*inv2 -
                               v*r*(-d2permw_dudv*permg - dpermw_du*dpermg_dv + dpermw_dv*dpermg_du + permw*d2permg_dudv)*inv2 +
                               2.0*v*r*(-dpermw_du*permg + permw*dpermg_du)*inv3*(dpermw_dv + dr_dv*permg + r*dpermg_dv)
                     );

                m.set(1, 1, 0, (2.0*permw*dpermw_du + dpermw_du*r*permg + permw*r*dpermg_du + 
                                v*d2permw_dudv*r*permg + v*dpermw_dv*r*dpermg_du -
                                v*dpermw_du*dr_dv*permg - v*permw*dr_dv*dpermg_du -
                                v*dpermw_du*r*dpermg_dv - v*permw*r*d2permg_dudv
                               )*inv2 -
                               2.0*(permw*permw + permw*r*permg + v*dpermw_dv*r*permg -
                                    v*permw*dr_dv*permg - v*permw*r*dpermg_dv
                                   )*inv3*(dpermw_du + r*dpermg_du)
                     );


                m.set(1, 1, 1, (2.0*permw*dpermw_dv + 2.0*dpermw_dv*r*permg + v*d2permw_dv2*r*permg -
                                v*permw*d2r_dv2*permg - 2.0*v*permw*dr_dv*dpermg_dv - v*permw*r*d2permg_dv2
                               )*inv2 -
                               2.0*(permw*permw + permw*r*permg + v*dpermw_dv*r*permg -
                                    v*permw*dr_dv*permg - v*permw*r*dpermg_dv
                                   )*inv3*(dpermw_dv + dr_dv*permg + r*dpermg_dv)
                     );
            }
        }
    }

    return 2;
}

IonFlux * IonFlux::clone() const {
    return new IonFlux(ratios, permeability);
}


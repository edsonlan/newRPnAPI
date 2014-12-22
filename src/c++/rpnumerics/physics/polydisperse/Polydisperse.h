#ifndef _POLYDISPERSE_
#define _POLYDISPERSE_

#include <stdio.h>
#include <stdlib.h>
#include "FluxFunction.h"

#include "Polydisperse_Params.h"

class Polydisperse : public FluxFunction {
    private:
//        double phimax;
//        double V1inf, V2inf;
//        double n1, n2;

        // Slip velocities
//        int Slip_jet(double phi1, double phi2, JetMatrix &uj, int degree) const;

        // Hindered settling factor
//        int Settling_jet(double phi1, double phi2, JetMatrix &vj, int degree) const;

    protected:
    public:
        Polydisperse(const Polydisperse_Params &);
        Polydisperse(const Polydisperse &);
        Polydisperse * clone() const;

        ~Polydisperse();

        int Hindered_jet(double phi1, double phi2, JetMatrix &Vj, int degree) const;
        int Relative_jet(double phi1, double phi2, JetMatrix &uj, int degree) const;
        int Absolute_jet(double phi1, double phi2, JetMatrix &vj, int degree) const;

        int jet(const WaveState &Phi, JetMatrix &flux, int degree) const;

};

#endif // _POLYDISPERSE_

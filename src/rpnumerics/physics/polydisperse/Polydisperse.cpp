#include "Polydisperse.h"

Polydisperse::Polydisperse(const Polydisperse_Params &param) : FluxFunction(param){
/** The DEFAULT parameters must setted differently for each model
 *
 *  For the MLB model, the initial data is
 *  // Maximum packing concentration, component(0)
 *     phimax = 1.0
 *  // Particle densities, component(1) and component(2)
 *     rho1 = 2.0
 *     rho2 = 3.0
 *  // The square diameters, component(3) and component(4)
 *     d1^2 = 1.0
 *     d2^2 = 0.25
 *  // Velocity exponents: there were used equally, component(5) and component(6)
 *     n1 = n2 = 5.0
 * 
 *  For the Bassoon restriction, the initial data is
 *  // Maximum packing concentration, component(0)
 *     phimax = 1.0
 *  // There is no particle densities, component(1) and component(2)
 *     rho1 = rho2 = 0.0
 *  // Terminal settling velocities, component(3) and component(4)
 *     -V1inf = - 1.0
 *     -Vinf2 = - 0.5
 *  // Velocity exponents, component(5) and component(6)
 *     n1 = 3.0
 *     n2 = 4.0
**/
}

Polydisperse * Polydisperse::clone() const {
    return new Polydisperse(*this);
}

Polydisperse::Polydisperse(const Polydisperse & copy): FluxFunction(copy.fluxParams()){
}

Polydisperse::~Polydisperse(){
}

/**
 * The polydisperse model given here is MLB (Masliyah-Lockett-Bassoon) with general power for velocities.
 * A <trick> in order to have the Bassoon restriction is possible by setting the particle densities as 0,
 * and treating the square of the particle densities as minus the terminal settling velocities (which may
 * be called "Stokes velocities").
 *
 * Thus the fluxes are the relative concentration multiplied by the absolute viscosity, i.e.,
 *      f_i(Phi) = phi_i * v_i(Phi)     for     i = 1, 2, ..., n
 *      and     Phi := (phi_1, phi_2, ..., phi_n).
 * Only the case n = 2 were implemented.
 *
 * The Absolute Velocity depends on the Relative Velocity which depends on the Hindered Settling Velocity 
 * and the Total Density. We have normalized the "Stokes constant" mu and the fluid density rho_f, so the
 * particle densities are actually normalized by the rho_f.
 *
 * Each one of these Velocities is implemented as a JET, so further modification would be easier.
 * (Pablo 28/02/2012)
**/

int Polydisperse::Hindered_jet(double phi1, double phi2, JetMatrix &Vj, int degree) const{
    double phi = phi1 + phi2;      // Total concentration

    double phimax = fluxParams().component(0);
    double n1 = fluxParams().component(5);
    double n2 = fluxParams().component(6);

    if (degree >= 0){
        // The velocity exponents must be larger than 1.0
        // It is not usual to have states outside the domain, however, we let this to happen.
        // The sign is controled inside the paramenter "inverse".
        //
        double modulus = fabs(1.0 - phi);
        double inverse = ( (phi != 1.0) ? 1./(1.0 - phi) : 0.0);
        double V1 = ( (phi < phimax) ? inverse*pow(modulus, n1 - 1.0) : 0.0 );
        double V2 = ( (phi < phimax) ? inverse*pow(modulus, n2 - 1.0) : 0.0 );

        Vj.set(0, V1);
        Vj.set(1, V2);

        if (degree >= 1){
            // We do not need the sign of the inverse, so we take the absolute value.
            inverse = fabs(inverse);

            double dV1_dphi = ( (phi < phimax) ? - (n1 - 2.0) * inverse * V1 : 0.0 );
            double dV2_dphi = ( (phi < phimax) ? - (n2 - 2.0) * inverse * V2 : 0.0 );

            Vj.set(0, 0, dV1_dphi);
//            Vj(0, 1, dV1_dphi);
            Vj.set(1, 0, dV2_dphi);
//            Vj(1, 1, dV2_dphi);

            if (degree >= 2){
                double d2V1_dphi2 = ( (phi < phimax) ? - (n1 - 3.0) * inverse * dV1_dphi : 0.0 );
                double d2V2_dphi2 = ( (phi < phimax) ? - (n2 - 3.0) * inverse * dV2_dphi : 0.0 );

                Vj.set(0, 0, 0, d2V1_dphi2);
//                Vj(0, 0, 1, d2V1_dphi2);
//                Vj(0, 1, 0, d2V1_dphi2);
//                Vj(0, 1, 1, d2V1_dphi2);
                Vj.set(1, 0, 0, d2V2_dphi2);
//                Vj(1, 0, 1, d2V2_dphi2);
//                Vj(1, 1, 0, d2V2_dphi2);
//                Vj(1, 1, 1, d2V2_dphi2);
            }
            else return -1; // ABORTED_PROCEDURE
        }
    }
    return 2; //SUCCESSFUL_PROCEDURE;
}

int Polydisperse::Relative_jet(double phi1, double phi2, JetMatrix &uj, int degree) const{
    double phi = phi1 + phi2;      // Total concentration

    double rho1 = fluxParams().component(1);
    double rho2 = fluxParams().component(2);
    double D1 = fluxParams().component(3);
    double D2 = fluxParams().component(4);

    JetMatrix hinderedj(2);
    Hindered_jet(phi1, phi2, hinderedj, degree);

    if (degree >= 0){
        double V1 = hinderedj.get(0);
        double V2 = hinderedj.get(1);

        // The TOTAL DENSITY is used along, recall that it is the normalized form
        double RHO = rho1*phi1 + rho2*phi2 + (1.0 - phi);

        double u1 = D1 * (rho1 - RHO) * V1;
        double u2 = D2 * (rho2 - RHO) * V2;

        uj.set(0, u1);
        uj.set(1, u2);

        if (degree >= 1){
            double dRHO_dphi1 = rho1 - 1.0;
            double dRHO_dphi2 = rho2 - 1.0;
            double dV1_dphi = hinderedj.get(0, 0);
            double dV2_dphi = hinderedj.get(1, 0);

            double du1_dphi1 = D1 * ( (rho1 - RHO) * dV1_dphi - dRHO_dphi1 * V1 );
            double du1_dphi2 = D1 * ( (rho1 - RHO) * dV1_dphi - dRHO_dphi2 * V1 );
            double du2_dphi1 = D2 * ( (rho2 - RHO) * dV2_dphi - dRHO_dphi1 * V2 );
            double du2_dphi2 = D2 * ( (rho2 - RHO) * dV2_dphi - dRHO_dphi2 * V2 );

            uj.set(0, 0, du1_dphi1);
            uj.set(0, 1, du1_dphi2);
            uj.set(1, 0, du2_dphi1);
            uj.set(1, 1, du2_dphi2);

            if (degree >= 2){
                double d2V1_dphi2 = hinderedj.get(0, 0, 0);
                double d2V2_dphi2 = hinderedj.get(1, 0, 0);

                double d2u1_dphi12 = D1 * ( (rho1 - RHO) * d2V1_dphi2 - 2.0 * dRHO_dphi1 * dV1_dphi );
                double d2u1_dphi22 = D1 * ( (rho1 - RHO) * d2V1_dphi2 - 2.0 * dRHO_dphi2 * dV1_dphi );

                double d2u2_dphi12 = D2 * ( (rho2 - RHO) * d2V2_dphi2 - 2.0 * dRHO_dphi1 * dV2_dphi );
                double d2u2_dphi22 = D2 * ( (rho2 - RHO) * d2V2_dphi2 - 2.0 * dRHO_dphi2 * dV2_dphi );

                double d2u1_dphi1phi2 = D1 * ( (rho1 - RHO) * d2V1_dphi2 + (2.0 - rho1 - rho2) * dV1_dphi );
                double d2u2_dphi1phi2 = D2 * ( (rho2 - RHO) * d2V2_dphi2 + (2.0 - rho1 - rho2) * dV2_dphi );

                uj.set(0, 0, 0, d2u1_dphi12   );
                uj.set(0, 1, 0, d2u1_dphi1phi2);
                uj.set(0, 0, 1, d2u1_dphi1phi2);
                uj.set(0, 1, 1, d2u1_dphi22   );
                uj.set(1, 0, 0, d2u2_dphi12   );
                uj.set(1, 1, 0, d2u2_dphi1phi2);
                uj.set(1, 0, 1, d2u2_dphi1phi2);
                uj.set(1, 1, 1, d2u2_dphi22   );
            }
            else return -1; // ABORTED_PROCEDURE
        }
    }
    return 2; //SUCCESSFUL_PROCEDURE;
}

int Polydisperse::Absolute_jet(double phi1, double phi2, JetMatrix &vj, int degree) const{
    JetMatrix relativej(2);
    Relative_jet(phi1, phi2, relativej, degree);

    if (degree >= 0){
        double u1 = relativej.get(0);
        double u2 = relativej.get(1);

        double v1 = (1.0 - phi1) * u1 - phi2 * u2;
        double v2 = (1.0 - phi2) * u2 - phi1 * u1;

        vj.set(0, v1);
        vj.set(1, v2);

        if (degree >= 1){
            double du1_dphi1 = relativej.get(0, 0);
            double du1_dphi2 = relativej.get(0, 1);
            double du2_dphi1 = relativej.get(1, 0);
            double du2_dphi2 = relativej.get(1, 1);

            double dv1_dphi1 = (1.0 - phi1) * du1_dphi1 - phi2 * du2_dphi1 - u1;
            double dv1_dphi2 = (1.0 - phi1) * du1_dphi2 - phi2 * du2_dphi2 - u2;

            double dv2_dphi1 = (1.0 - phi2) * du2_dphi1 - phi1 * du1_dphi1 - u1;
            double dv2_dphi2 = (1.0 - phi2) * du2_dphi2 - phi1 * du1_dphi2 - u2;

            vj.set(0, 0, dv1_dphi1);
            vj.set(0, 1, dv1_dphi2);

            vj.set(1, 0, dv2_dphi1);
            vj.set(1, 1, dv2_dphi2);

            if (degree >= 2){
                double d2u1_dphi12    = relativej.get(0, 0, 0);
                double d2u1_dphi1phi2 = relativej.get(0, 0, 1);
                double d2u1_dphi22    = relativej.get(0, 1, 1);
                double d2u2_dphi12    = relativej.get(1, 0, 0);
                double d2u2_dphi1phi2 = relativej.get(1, 0, 1);
                double d2u2_dphi22    = relativej.get(1, 1, 1);

                double d2v1_dphi12    = (1.0 - phi1) * d2u1_dphi12 - phi2 * d2u2_dphi12 - 2.0 * du1_dphi1;
                double d2v1_dphi1phi2 = (1.0 - phi1) * d2u1_dphi1phi2 - phi2 * d2u2_dphi1phi2 - du1_dphi2 - du2_dphi1;
                double d2v1_dphi22    = (1.0 - phi1) * d2u1_dphi22 - phi2 * d2u2_dphi22 - 2.0 * du2_dphi2;

                double d2v2_dphi12    = (1.0 - phi2) * d2u2_dphi12 - phi1 * d2u1_dphi12 - 2.0 * du1_dphi1;
                double d2v2_dphi1phi2 = (1.0 - phi2) * d2u2_dphi1phi2 - phi1 * d2u1_dphi1phi2 - du1_dphi2 - du2_dphi1;
                double d2v2_dphi22    = (1.0 - phi2) * d2u2_dphi22 - phi1 * d2u1_dphi22 - 2.0 * du2_dphi2;

                vj.set(0, 0, 0, d2v1_dphi12   );
                vj.set(0, 0, 1, d2v1_dphi1phi2);
                vj.set(0, 1, 0, d2v1_dphi1phi2);
                vj.set(0, 1, 1, d2v1_dphi22   );

                vj.set(1, 0, 0, d2v2_dphi12   );
                vj.set(1, 0, 1, d2v2_dphi1phi2);
                vj.set(1, 1, 0, d2v2_dphi1phi2);
                vj.set(1, 1, 1, d2v2_dphi22   );
            }
            else return -1; // ABORTED_PROCEDURE
        }
    }
    return 2; //SUCCESSFUL_PROCEDURE;
}

int Polydisperse::jet(const WaveState &Phi, JetMatrix &flux, int degree) const{
    double phi1 = Phi(0);
    double phi2 = Phi(1);

    JetMatrix absolutej(2);
    Absolute_jet(phi1, phi2, absolutej, degree);

    if (degree >= 0){
        double v1 = absolutej.get(0);
        double v2 = absolutej.get(1);
        double f1 = phi1 * v1;
        double f2 = phi2 * v2;

        flux.set(0, f1);
        flux.set(1, f2);

        if (degree >= 1){
            double dv1_dphi1 = absolutej.get(0, 0);
            double dv1_dphi2 = absolutej.get(0, 1);
            double dv2_dphi1 = absolutej.get(1, 0);
            double dv2_dphi2 = absolutej.get(1, 1);

            double df1_dphi1 = v1 + phi1 * dv1_dphi1;
            double df1_dphi2 = phi1 * dv1_dphi2;

            double df2_dphi1 = phi2 * dv2_dphi1;
            double df2_dphi2 = v2 + phi2 * dv2_dphi2;

            flux.set(0, 0, df1_dphi1);
            flux.set(0, 1, df1_dphi2);
            flux.set(1, 0, df2_dphi1);
            flux.set(1, 1, df2_dphi2);

            if (degree >= 2){
                double d2v1_dphi12    = absolutej.get(0, 0, 0);
                double d2v1_dphi1phi2 = absolutej.get(0, 0, 1);
                double d2v1_dphi22    = absolutej.get(0, 1, 1);
                double d2v2_dphi12    = absolutej.get(1, 0, 0);
                double d2v2_dphi1phi2 = absolutej.get(1, 0, 1);
                double d2v2_dphi22    = absolutej.get(1, 1, 1);

                double d2f1_dphi12    = 2.0 * dv1_dphi1 + phi1 * d2v1_dphi12;
                double d2f1_dphi1phi2 = dv1_dphi2 + phi1 * d2v1_dphi1phi2;
                double d2f1_dphi22    = phi1 * d2v1_dphi22;

                double d2f2_dphi12    = phi2 * d2v2_dphi12;
                double d2f2_dphi1phi2 = dv2_dphi1 + phi2 * d2v2_dphi1phi2;
                double d2f2_dphi22    = 2.0 * dv2_dphi2 + phi2 * d2v2_dphi22;

                flux.set(0, 0, 0, d2f1_dphi12   );
                flux.set(0, 0, 1, d2f1_dphi1phi2);
                flux.set(0, 1, 0, d2f1_dphi1phi2);
                flux.set(0, 1, 1, d2f1_dphi22   );
                flux.set(1, 0, 0, d2f2_dphi12   );
                flux.set(1, 0, 1, d2f2_dphi1phi2);
                flux.set(1, 1, 0, d2f2_dphi1phi2);
                flux.set(1, 1, 1, d2f2_dphi22   );
            }
            else return -1; // ABORTED_PROCEDURE
        }
    }
    return 2; //SUCCESSFUL_PROCEDURE;
}

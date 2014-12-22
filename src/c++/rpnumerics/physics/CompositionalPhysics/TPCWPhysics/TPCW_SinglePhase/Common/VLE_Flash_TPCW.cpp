#include "VLE_Flash_TPCW.h"

VLE_Flash_TPCW::VLE_Flash_TPCW(const MolarDensity *mdl_, const MolarDensity *mdv_){
    mdl = (MolarDensity*)mdl_;
    mdv = (MolarDensity*)mdv_;
    
    
    std::cout<<"mdl"<<mdl<<"  "<<mdv<<std::endl;

    Epsilon_machine = pow(0.1, 10); // Era 10^(-14)

    // Initialization of the molar fractions : this may change for the complete FLASH method
    zc = 0.5;
    zw = 0.5;

    R  = mdl->get_R();
    P  = mdl->get_P();
    q1 = mdl->get_q1();
    q2 = mdl->get_q2();

    bc = mdl->get_bc();
    bw = mdl->get_bw();

    Pcc = mdl->get_Pcc();
    Pcw = mdl->get_Pcw();
    Tcc = mdl->get_Tcc();
    Tcw = mdl->get_Tcw();

    omega_c = mdl->get_omega_c();
    omega_w = mdl->get_omega_w();

    ROOTTWO = 1.414213562373095; // Square root of two
    INV_2R2 = 0.353553390593274; // Inverse of 2 times square root of two
}

VLE_Flash_TPCW::~VLE_Flash_TPCW(){}

// Fugacity Liquid CO2
//
int VLE_Flash_TPCW::fugacitycoef_lc(double x, double T, int degree, JetMatrix &fugacitycoef_lcj){
    double bi = bc;

    JetMatrix bj(1);    //printf("fugacitycoef_lc1: mdl = %p\n", mdl);

    mdl->b_jet(x, degree, bj);


    JetMatrix epsilonj(2);
    mdl->Epsilon_jet(x, T, degree, epsilonj);

    JetMatrix epsilon_ij(1);
    mdl->epsilon_c_jet(T, degree, epsilon_ij);

    JetMatrix gammaj(2);
    mdl->gamma1_jet(x, T, degree, gammaj);

    // Construction of epsilon_bar
    JetMatrix epsilon_barj(2);
    epsilon_bar_jet(x, T, bi, degree, bj, epsilonj, epsilon_ij, gammaj, epsilon_barj);

    JetMatrix Bj(2);
    mdl->B_jet(x, T, degree, Bj);

    JetMatrix Zj(2);
    mdl->Z_jet(x, T, degree, Zj);

    return fugacity_coeficient_jet(x, T, bi, degree, bj, Bj, Zj, epsilon_barj, fugacitycoef_lcj);
}

// Fugacity Vapor CO2
//
int VLE_Flash_TPCW::fugacitycoef_vc(double x, double T, int degree, JetMatrix &fugacitycoef_vcj){
    double bi = bc;

    JetMatrix bj(1);
    mdv->b_jet(x, degree, bj);

    JetMatrix epsilonj(2);
    mdv->Epsilon_jet(x, T, degree, epsilonj);

    JetMatrix epsilon_ij(1);
    mdv->epsilon_c_jet(T, degree, epsilon_ij);

    JetMatrix gammaj(2);
    mdv->gamma1_jet(x, T, degree, gammaj);

    // Construction of epsilon_bar
    JetMatrix epsilon_barj(2);
    epsilon_bar_jet(x, T, bi, degree, bj, epsilonj, epsilon_ij, gammaj, epsilon_barj);

    JetMatrix Bj(2);
    mdv->B_jet(x, T, degree, Bj);

    JetMatrix Zj(2);
    mdv->Z_jet(x, T, degree, Zj);

    return fugacity_coeficient_jet(x, T, bi, degree, bj, Bj, Zj, epsilon_barj, fugacitycoef_vcj);
}

// Fugacity Liquid H2O
//
int VLE_Flash_TPCW::fugacitycoef_lw(double x, double T, int degree, JetMatrix &fugacitycoef_lwj){
    double bi = bw;

    JetMatrix bj(1);
    mdl->b_jet(x, degree, bj);

    JetMatrix epsilonj(2);
    mdl->Epsilon_jet(x, T, degree, epsilonj);

    JetMatrix epsilon_ij(1);
    mdl->epsilon_w_jet(T, degree, epsilon_ij);

    JetMatrix gammaj(2);
    mdl->gamma2_jet(x, T, degree, gammaj);

    // Construction of epsilon_bar
    JetMatrix epsilon_barj(2);
    epsilon_bar_jet(x, T, bi, degree, bj, epsilonj, epsilon_ij, gammaj, epsilon_barj);

    JetMatrix Bj(2);
    mdl->B_jet(x, T, degree, Bj);

    JetMatrix Zj(2);
    mdl->Z_jet(x, T, degree, Zj);

    return fugacity_coeficient_jet(x, T, bi, degree, bj, Bj, Zj, epsilon_barj, fugacitycoef_lwj);
}

// Fugacity Vapor H2O
//
int VLE_Flash_TPCW::fugacitycoef_vw(double x, double T, int degree, JetMatrix &fugacitycoef_vwj){
    double bi = bw;

    JetMatrix bj(1);
    mdv->b_jet(x, degree, bj);

    JetMatrix epsilonj(2);
    mdv->Epsilon_jet(x, T, degree, epsilonj);

    JetMatrix epsilon_ij(1);
    mdv->epsilon_w_jet(T, degree, epsilon_ij);

    JetMatrix gammaj(2);
    mdv->gamma2_jet(x, T, degree, gammaj);

    // Construction of epsilon_bar
    JetMatrix epsilon_barj(2);
    epsilon_bar_jet(x, T, bi, degree, bj, epsilonj, epsilon_ij, gammaj, epsilon_barj);

    JetMatrix Bj(2);
    mdv->B_jet(x, T, degree, Bj);

    JetMatrix Zj(2);
    mdv->Z_jet(x, T, degree, Zj);

    return fugacity_coeficient_jet(x, T, bi, degree, bj, Bj, Zj, epsilon_barj, fugacitycoef_vwj);
}


// Fugacity coeficient for different types
//
int VLE_Flash_TPCW::fugacity_coeficient_jet(double x, double T, double bi, int degree,
                    JetMatrix &bj, JetMatrix &Bj, JetMatrix &Zj, JetMatrix &epsilon_barj,
                    JetMatrix &fugacity_coeficientj){

    if (degree >= 0){
        double b = bj.get(0);
        double B = Bj.get(0);
        double Z = Zj.get(0); // Compressibility Factor

        double E = epsilon_barj.get(0);
        double inv_ZB = 1.0/(Z - B);
        double log_ZB = log((Z + (1.0 + ROOTTWO)*B)/(Z + (1.0 - ROOTTWO)*B));

        // This is the logarithm of the fugacity coeficient [ + log(Z - B) ]:
        double log_fc = bi*(Z - 1.0)/b - INV_2R2*E*log_ZB;
        double fc = exp(log_fc)*inv_ZB;

        fugacity_coeficientj.set(0, fc);

        if (degree >= 1){
            double db_dx = bj.get(0,0);
            double dB_dx = Bj.get(0,0);
            double dB_dT = Bj.get(0,1);
            double dZ_dx = Zj.get(0,0); // Compressibility Factor
            double dZ_dT = Zj.get(0,1);

            double dE_dx = epsilon_barj.get(0,0);
            double dE_dT = epsilon_barj.get(0,1);

            double base_ZB = 1.0/( (Z + (1.0 + ROOTTWO)*B)*(Z + (1.0 - ROOTTWO)*B) );

            double dlogfc_dx = bi*(dZ_dx - db_dx*(Z - 1.0)/b)/b - INV_2R2*dE_dx*log_ZB + E*(B*dZ_dx - Z*dB_dx)*base_ZB;
            double dlogfc_dT = bi*dZ_dT/b                       - INV_2R2*dE_dT*log_ZB + E*(B*dZ_dT - Z*dB_dT)*base_ZB;

            double interior_x = (dB_dx - dZ_dx)*inv_ZB + dlogfc_dx;
            double interior_T = (dB_dT - dZ_dT)*inv_ZB + dlogfc_dT;

            double dfc_dx = interior_x*fc;
            double dfc_dT = interior_T*fc;

            fugacity_coeficientj.set(0, 0, dfc_dx);
            fugacity_coeficientj.set(0, 1, dfc_dT);

            if (degree == 2){ //Testar segundas derivadas respeito de T.

                double d2b_dx2  = bj.get(0,0,0);
                double d2B_dx2  = Bj.get(0,0,0);
                double d2B_dxdT = Bj.get(0,0,1);
                double d2B_dT2  = Bj.get(0,1,1);
                double d2Z_dx2  = Zj.get(0,0,0); // Compressibility Factor
                double d2Z_dxdT = Zj.get(0,0,1);
                double d2Z_dT2  = Zj.get(0,1,1);

                double d2E_dx2  = epsilon_barj.get(0,0,0);
                double d2E_dxdT = epsilon_barj.get(0,0,1);
                double d2E_dT2  = epsilon_barj.get(0,1,1);

                double inv_ZB2  = inv_ZB * inv_ZB;
                double base_ZB2 = base_ZB*base_ZB;

                double d2logfc_dx2  = bi*(d2Z_dx2*b - d2b_dx2*(Z - 1.0) - 2.0*dZ_dx*db_dx)/(b*b) + 2.0*bi*db_dx*db_dx*(Z - 1.0)/(b*b*b)
                                    - INV_2R2*d2E_dx2*log_ZB + 2.0*dE_dx*(B*dZ_dx - Z*dB_dx)*base_ZB
                                    + E*(B*d2Z_dx2 - Z*d2B_dx2)*base_ZB
                                    - 2.0*E*(B*dZ_dx - Z*dB_dx)*(Z*dZ_dx + Z*dB_dx + B*dZ_dx - B*dB_dx)*base_ZB2;
                double d2logfc_dxdT = bi*(d2Z_dxdT*b - db_dx*dZ_dT)/(b*b) - INV_2R2*d2E_dxdT*log_ZB
                                    + ( dE_dT*(B*dZ_dx - Z*dB_dx) + dE_dx*(B*dZ_dT - Z*dB_dT) )*base_ZB
                                    + E*(dB_dx*dZ_dT + B*d2Z_dxdT - Z*d2B_dxdT - dB_dT*dZ_dx)*base_ZB
                                    - 2.0*E*(B*dZ_dT - Z*dB_dT)*(Z*dZ_dx + Z*dB_dx + B*dZ_dx - B*dB_dx)*base_ZB2;
                double d2logfc_dT2  = bi*d2Z_dT2/b - INV_2R2*d2E_dT2*log_ZB + 2.0*dE_dT*(B*dZ_dT - Z*dB_dT)*base_ZB
                                    + E*(B*d2Z_dT2 - Z*d2B_dT2)*base_ZB
                                    - 2.0*E*(B*dZ_dT - Z*dB_dT)*(Z*dZ_dT + Z*dB_dT + B*dZ_dT - B*dB_dT)*base_ZB2;

                double d2fc_dx2  = ( (d2B_dx2  -  d2Z_dx2)*inv_ZB +(dB_dx - dZ_dx)*(dB_dx - dZ_dx)*inv_ZB2 + d2logfc_dx2  + interior_x*interior_x )*fc;
                double d2fc_dxdT = ( (d2B_dxdT - d2Z_dxdT)*inv_ZB +(dB_dx - dZ_dx)*(dB_dT - dZ_dT)*inv_ZB2 + d2logfc_dxdT + interior_x*interior_T )*fc;
                double d2fc_dTdx = d2fc_dxdT;
                double d2fc_dT2  = ( (d2B_dT2  -  d2Z_dT2)*inv_ZB +(dB_dT - dZ_dT)*(dB_dT - dZ_dT)*inv_ZB2 + d2logfc_dT2  + interior_T*interior_T )*fc;

                fugacity_coeficientj.set(0, 0, 0, d2fc_dx2);
                fugacity_coeficientj.set(0, 0, 1, d2fc_dxdT);
                fugacity_coeficientj.set(0, 1, 0, d2fc_dTdx);
                fugacity_coeficientj.set(0, 1, 1, d2fc_dT2);

            }
            else return -1; // ABORTED_PROCEDURE
        }
    }

    return 2; // SUCCESSFUL_PROCEDURE
}


int VLE_Flash_TPCW::epsilon_bar_jet(double x, double T, double bi, int degree,
                    JetMatrix &bj, JetMatrix &epsilonj, JetMatrix &epsilon_ij, JetMatrix &gammaj,
                    JetMatrix &epsilon_barj){

    if (degree >= 0){
        double gamma = gammaj.get(0);
        double eps   = epsilonj.get(0);
        double eps_i = epsilon_ij.get(0);
        double b     = bj.get(0);

        // Auxiliary expressions for the JET
        double b_inv    = 1.0/b;
        double first_part = q1*eps_i + q2*(eps*eps + eps_i*eps_i) + log(gamma) + log(b/bi) + bi*b_inv - 1.0;
        double qs_inv     = 1.0/(q1 + 2.0*q2*eps);

        double epsilon_bar = first_part*qs_inv;

        epsilon_barj.set(0, epsilon_bar);

        if (degree >= 1){
            double dgamma_dx = gammaj.get(0,0);
            double dgamma_dT = gammaj.get(0,1);

            double deps_dx = epsilonj.get(0,0);
            double deps_dT = epsilonj.get(0,1);

            double depsi_dT = epsilon_ij.get(0,0);
            double db_dx    = bj.get(0,0);

            // More auxiliary expressions for the JET
            double b2_inv    = b_inv * b_inv;
            double qs2_inv   = qs_inv*qs_inv;
            double gamma_inv = 1.0/gamma;

            // More auxiliary derivative expressions for the JET
            double dfirst_dx = 2.0*q2*eps*deps_dx + gamma_inv*dgamma_dx + (b_inv - bi*b2_inv)*db_dx;
            double dfirst_dT = (q1 + 2.0*q2*eps_i)*depsi_dT + 2.0*q2*eps*deps_dT + gamma_inv*dgamma_dT;

            double depsilon_bar_dx = dfirst_dx*qs_inv - 2.0*q2*qs2_inv*deps_dx*first_part;
            double depsilon_bar_dT = dfirst_dT*qs_inv - 2.0*q2*qs2_inv*deps_dT*first_part;

            epsilon_barj.set(0, 0, depsilon_bar_dx);
            epsilon_barj.set(0, 1, depsilon_bar_dT);

            if (degree == 2){
                double d2gamma_dx2  = gammaj.get(0,0,0);
                double d2gamma_dxdT = gammaj.get(0,0,1);
                double d2gamma_dTdx = gammaj.get(0,1,0);
                double d2gamma_dT2  = gammaj.get(0,1,1);

                double d2eps_dx2  = epsilonj.get(0,0,0);
                double d2eps_dxdT = epsilonj.get(0,0,1);
                double d2eps_dTdx = epsilonj.get(0,1,0);
                double d2eps_dT2  = epsilonj.get(0,1,1);

                double d2epsi_dT2 = epsilon_ij.get(0,0,0);
                double d2b_dx2    = bj.get(0,0,0);

                // More auxiliary expressions for the JET
                double b3_inv     = b2_inv * b_inv;
                double qs3_inv    = qs2_inv*qs_inv;
                double gamma2_inv = gamma_inv*gamma_inv;

                // More auxiliary derivative expressions for the JET
                double d2first_dx2  = 2.0*q2*eps*d2eps_dx2  + 2.0*q2*deps_dx*deps_dx - gamma2_inv*dgamma_dx*dgamma_dx
                                    + gamma_inv*d2gamma_dx2 - (b2_inv - 2.0*bi*b3_inv)*db_dx*db_dx + (b_inv - bi*b2_inv)*d2b_dx2;
                double d2first_dxdT = 2.0*q2*eps*d2eps_dxdT + 2.0*q2*deps_dx*deps_dT - gamma2_inv*dgamma_dx*dgamma_dT
                                    + gamma_inv*d2gamma_dxdT;
                double d2first_dT2  = (q1 + 2.0*q2*eps_i)*d2epsi_dT2  + 2.0*q2*depsi_dT*depsi_dT - gamma2_inv*dgamma_dT*dgamma_dT
                                    + gamma_inv*d2gamma_dT2 + 2.0*q2*deps_dT*deps_dT + 2.0*q2*eps*d2eps_dT2;

                double d2qs_inv_dx2  = 8.0*q2*q2*qs3_inv*deps_dx*deps_dx - 2.0*q2*qs2_inv*d2eps_dx2;
                double d2qs_inv_dxdT = 8.0*q2*q2*qs3_inv*deps_dx*deps_dT - 2.0*q2*qs2_inv*d2eps_dxdT;
                double d2qs_inv_dT2  = 8.0*q2*q2*qs3_inv*deps_dT*deps_dT - 2.0*q2*qs2_inv*d2eps_dT2;

                double d2epsilon_bar_dx2  = d2first_dx2 *qs_inv + first_part*d2qs_inv_dx2  - 4.0*q2*qs2_inv*deps_dx*dfirst_dx;
                double d2epsilon_bar_dxdT = d2first_dxdT*qs_inv + first_part*d2qs_inv_dxdT - 2.0*q2*qs2_inv*(deps_dx*dfirst_dT
                                          + deps_dT*dfirst_dx);
                double d2epsilon_bar_dTdx = d2epsilon_bar_dxdT;
                double d2epsilon_bar_dT2  = d2first_dT2 *qs_inv + first_part*d2qs_inv_dT2  - 4.0*q2*qs2_inv*deps_dT*dfirst_dT;

                epsilon_barj.set(0, 0, 0, d2epsilon_bar_dx2);
                epsilon_barj.set(0, 0, 1, d2epsilon_bar_dxdT);
                epsilon_barj.set(0, 1, 0, d2epsilon_bar_dTdx);
                epsilon_barj.set(0, 1, 1, d2epsilon_bar_dT2);
            }
            else return -1; // ABORTED_PROCEDURE
        }
    }

    return 2; // SUCCESSFUL_PROCEDURE
}

int VLE_Flash_TPCW::molarfractions_jet(double T, JetMatrix &xcj, JetMatrix &ywj, int degree){

    if (degree >= 0){
        double xc, yw;

        flash(T, xc, yw);
        xcj.set(0, xc);
        ywj.set(0, yw);

        if (degree >= 1){
            JetMatrix fugacitycoef_lcj(2);
            JetMatrix fugacitycoef_vcj(2);
            JetMatrix fugacitycoef_lwj(2);
            JetMatrix fugacitycoef_vwj(2);

            fugacitycoef_lc(xc, T, degree, fugacitycoef_lcj);
            fugacitycoef_vc(yw, T, degree, fugacitycoef_vcj);
            fugacitycoef_lw(xc, T, degree, fugacitycoef_lwj);
            fugacitycoef_vw(yw, T, degree, fugacitycoef_vwj);

            double phicL = fugacitycoef_lcj.get(0);
            double phicV = fugacitycoef_vcj.get(0);
            double phiwL = fugacitycoef_lwj.get(0);
            double phiwV = fugacitycoef_vwj.get(0);

            double dphicL_dx = fugacitycoef_lcj.get(0,0);
            double dphicV_dx = fugacitycoef_vcj.get(0,0);
            double dphiwL_dx = fugacitycoef_lwj.get(0,0);
            double dphiwV_dx = fugacitycoef_vwj.get(0,0);

            double dphicL_dT = fugacitycoef_lcj.get(0,1);
            double dphicV_dT = fugacitycoef_vcj.get(0,1);
            double dphiwL_dT = fugacitycoef_lwj.get(0,1);
            double dphiwV_dT = fugacitycoef_vwj.get(0,1);

            // We solve the linear problem Ax = b in order to find the molar fraction derivatives.
            double A11 = phicL + xc*dphicL_dx;
            double A12 = phicV - (1.0 - yw)*dphicV_dx;
            double A21 = phiwL - (1.0 - xc)*dphiwL_dx;
            double A22 = phiwV + yw*dphiwV_dx;

            double b1 = (1.0 - yw)*dphicV_dT - xc*dphicL_dT;
            double b2 = (1.0 - xc)*dphiwL_dT - yw*dphiwV_dT;

            double discriminant = A11*A22 - A21*A12;
            double Discr_mean_2 = (A11 + A12 + A21 + A22)*(A11 + A12 + A21 + A22)/16.0;
            if ( fabs(discriminant) < Epsilon_machine*Discr_mean_2 ){
                printf("The VLE matrix is near singular with xc = %lf, yw = %lf at temperature %.4lf K degrees\n", xc, yw, T);
                return -1; // ABORTED_PROCEDURE
            }
            double inv_discr = 1.0/discriminant;

            double dxc_dT = (b1*A22 - b2*A12) * inv_discr;
            double dyw_dT = (A11*b2 - A21*b1) * inv_discr;

            xcj.set(0, 0, dxc_dT);
            ywj.set(0, 0, dyw_dT);

            if (degree == 2){
                double d2phicL_dx2  = fugacitycoef_lcj.get(0,0,0);
                double d2phicV_dx2  = fugacitycoef_vcj.get(0,0,0);
                double d2phiwL_dx2  = fugacitycoef_lwj.get(0,0,0);
                double d2phiwV_dx2  = fugacitycoef_vwj.get(0,0,0);

                double d2phicL_dxdT = fugacitycoef_lcj.get(0,1,0);
                double d2phicV_dxdT = fugacitycoef_vcj.get(0,1,0);
                double d2phiwL_dxdT = fugacitycoef_lwj.get(0,1,0);
                double d2phiwV_dxdT = fugacitycoef_vwj.get(0,1,0);

                double d2phicL_dT2  = fugacitycoef_lcj.get(0,1,1);
                double d2phicV_dT2  = fugacitycoef_vcj.get(0,1,1);
                double d2phiwL_dT2  = fugacitycoef_lwj.get(0,1,1);
                double d2phiwV_dT2  = fugacitycoef_vwj.get(0,1,1);

                // Notice that the linear problema now is Ax = c, with the same matrix A.
                double c1 = (1.0 - yw)*(d2phicV_dT2 + 2.0*d2phicV_dxdT*dyw_dT + d2phicV_dx2*dyw_dT*dyw_dT)
                          - xc*(d2phicL_dT2 + 2.0*d2phicL_dxdT*dxc_dT + d2phicL_dx2*dxc_dT*dxc_dT)
                          - 2.0*dxc_dT*(dphicL_dT + dphicL_dx*dxc_dT) - 2.0*dyw_dT*(dphicV_dT + dphicV_dx*dyw_dT);
                double c2 = (1.0 - xc)*(d2phiwL_dT2 + 2.0*d2phiwL_dxdT*dxc_dT + d2phiwL_dx2*dxc_dT*dxc_dT)
                          - yw*(d2phiwV_dT2 + 2.0*d2phiwV_dxdT*dyw_dT + d2phiwV_dx2*dyw_dT*dyw_dT)
                          - 2.0*dxc_dT*(dphiwL_dT + dphiwL_dx*dxc_dT) - 2.0*dyw_dT*(dphiwV_dT + dphiwV_dx*dyw_dT);

                double d2xc_dT2 = (c1*A22 - c2*A12) * inv_discr;
                double d2yw_dT2 = (A11*c2 - A21*c1) * inv_discr;

                xcj.set(0, 0, 0, d2xc_dT2);
                ywj.set(0, 0, 0, d2yw_dT2);

            }
            else return -1; // ABORTED_PROCEDURE
        }
    }

    return 2; // SUCCESSFUL_PROCEDURE
}




int VLE_Flash_TPCW::flash(double T, double &xc, double &yw){
    //printf("flash1\n");

    double vapor_frac = 0.05;

    // First estimatives using the Wilson Equation, [Orr, Eq. (3.4.2)]
    // Estimates are correct, in agreement with Ehsan's correction suggested on the 19 Apr, 2011.
    
    // Remember Ki = yi/xi , for the components i=c,w.
    double Kc, Kw;

//    Kc = Pcc/P*exp(5.37*(1.0 + omega_c)*(1.0 - Tcc/T)); HERE WE CAN PUT A FLAG FOR CHOOSING EITHER WILSON OR POLYS ESTIMATES.
//    Kw = Pcw/P*exp(5.37*(1.0 + omega_w)*(1.0 - Tcw/T));

    // Here we define the new estimatives for Kc and Kw.  
    // These estimatives have found by previously running the FLASH calculation and by comparing
    // simple formula-interpolations using the CFTOOL by MATLAB.

    // For Kc we write Kc = ((1-yw)=yc)/xc WE INDICATE BELOW THAT THE VALUES ARE ESTIMATED FOR THE FIRST KICK.
    // For Kw we write Kw = yw/((1-xc)=xw) 

//    MATLAB OUTPUT FOR FITTINGS 
      double T2 = T*T;
      double T3 = T2*T;
      double T4 = T3*T;


//    Model Poly3:
//    Coefficients (with 95% confidence bounds):
      double a1 = -8.685e-09;       // (-8.767e-09, -8.602e-09)
      double a2 =  1.113e-05;       // (1.104e-05, 1.122e-05)
      double a3 = -0.004798;        // (-0.004831, -0.004764)
      double a4 =  0.703;           // (0.6989, 0.7071)
      double xc_est = a1*T3 + a2*T2 + a3*T + a4;

//    Model Poly4:
//    Coefficients (with 95% confidence bounds):
      double b1 =  2.362e-10;       // (2.313e-10, 2.41e-10)
      double b2 = -3.021e-07;       // (-3.093e-07, -2.949e-07)
      double b3 =  0.0001475;       // (0.0001435, 0.0001515)
      double b4 = -0.0325;          // (-0.03347, -0.03153)
      double b5 =  2.723;           // (2.635, 2.81)
      double yw_est = b1*T4 + b2*T3 + b3*T2 + b4*T + b5;


//  This is the first kick for the values of the equilibrium K-values.
    Kc = (1.0 - yw_est )/xc_est    ;
    Kw = yw_est / ( 1.0 - xc_est ) ;

    double F, dF_dv;
    double Err1, Err2, Err3;
    Err1 = Err2 = Err3 = 1.0;

    double xw, yc, sum_x, sum_y;
    double vapor_fugacity;
    JetMatrix liquid_fugj(2), vapor_fugj(2);

    int counter = 0; //printf("flash2\n");

    while( (Err1 > Epsilon_machine) || (Err2 > Epsilon_machine) || (Err3 > Epsilon_machine) ){
        F  = zc*(Kc - 1.0)/(1.0 + vapor_frac*(Kc - 1.0)) + zw*(Kw - 1.0)/(1.0 + vapor_frac*(Kw - 1.0));
        dF_dv = - ( zc*(Kc - 1.0)/(1.0 + vapor_frac*(Kc - 1.0))*(Kc - 1.0)/(1.0 + vapor_frac*(Kc - 1.0))
              + zw*(Kw - 1.0)/(1.0 + vapor_frac*(Kw - 1.0))*(Kw - 1.0)/(1.0 + vapor_frac*(Kw - 1.0)) );

        

        if(dF_dv == 0.0){
            return -1; // ABORTED_PROCEDURE
        }
        vapor_frac = vapor_frac - F/dF_dv;
        vapor_frac = clamp(vapor_frac, 0.0, 1.0);

        xc = zc/(1.0 + vapor_frac*(Kc - 1.0));
        xw = zw/(1.0 + vapor_frac*(Kw - 1.0));

        yc = Kc*xc;
        yw = Kw*xw;

        sum_x = xc + xw;
        sum_y = yc + yw;

        Err1 = fabs(sum_x - 1.0);
        Err2 = fabs(sum_y - 1.0);
        Err3 = F;

        // Normalization of the compositions
        xc = xc/sum_x;
        yw = yw/sum_y;

        

        // If we give an initial guess for xc and yw, these two procedures must be at the begining of the while loop
        fugacitycoef_lc(xc, T, 0, liquid_fugj);//printf("counter = %d; F = %g; DF = %g\n", counter, F, dF_dv);
        fugacitycoef_vc(yw, T, 0, vapor_fugj);
        vapor_fugacity = vapor_fugj.get(0);
        Kc = ( vapor_fugacity == 0.0 ? 0.0 : liquid_fugj.get(0)/vapor_fugacity );

        fugacitycoef_lw(xc, T, 0, liquid_fugj);
        fugacitycoef_vw(yw, T, 0, vapor_fugj);
        vapor_fugacity = vapor_fugj.get(0);
        Kw = ( vapor_fugacity == 0.0 ? 0.0 : liquid_fugj.get(0)/vapor_fugacity );

        if(counter > 50) {
            // The following prints may be removed later
            printf("Newton method exception (break) after %d iterations at %.4lf K degrees\n", counter, T);
            printf("Error1 = %20.20lf, ", Err1);
            printf("Error2 = %20.20lf, ", Err2);
            printf("Error3 = %20.20lf\n", Err3);
            return -2; // NEWTON_DOES_NOT_CONVERGE
        }
        counter++; // This counter avoids infinite loop
    }

    return 2; // SUCCESSFUL_PROCEDURE
}


#include "Thermodynamics.h"

double Thermodynamics::Tref_rock  = 273.15;
double Thermodynamics::Tref_water = 274.3775;
double Thermodynamics::P          = 100.9e5;

double Thermodynamics::Rock_Cr = 2.029e6;
double Thermodynamics::Cw_     = 4297.0;

// Typical values for the adimensionalization
double Thermodynamics::T_typical_   = 304.63;
double Thermodynamics::Rho_typical_ = 998.2;
double Thermodynamics::U_typical_   = 4.22e-3;
// h_typical_ = Cw_ * (T_typical_ - Tref_water);
double Thermodynamics::h_typical_ = 129994.9925;
double Thermodynamics::rhoW_const = 998.2;

// Some constants for mug
double Thermodynamics::a0 = -1.94760101098783e-6;
double Thermodynamics::a1 =  0.013524080086578;
double Thermodynamics::a2 = -9.043578102452411;
double Thermodynamics::a3 =  1.612763701298628e3;

// some constants for muw
double Thermodynamics::b0 = -0.0123274;
double Thermodynamics::b1 = 27.1038;
double Thermodynamics::b2 = -23527.5;
double Thermodynamics::b3 =  1.01425e7;
double Thermodynamics::b4 = -2.17342e9;
double Thermodynamics::b5 =  1.86935e11;

VLE_Flash_TPCW *Thermodynamics::Flash = (VLE_Flash_TPCW *)0;

// Store the compositions for a fixed temperature
double Thermodynamics::xc_T, Thermodynamics::yw_T;
double Thermodynamics::dxc_dTheta, Thermodynamics::d2xc_dTheta2;
double Thermodynamics::dyw_dTheta, Thermodynamics::d2yw_dTheta2;

double Thermodynamics::used_Theta;

using namespace std;
void Thermodynamics::SetCompositions(double Theta) {
    //std::cout << "Thermodynamics::SetCompositions: Theta = " << Theta << ", used_Theta = " << used_Theta << std::endl;

    if (used_Theta != Theta){
        //std::cout << "Thermodynamics::SetCompositions: new Theta. Flash = " << Flash << std::endl;
//        Flash->flash(Theta2T(Theta), xc_T, yw_T);

        // Flash tiene que rellenar esto de aqui abajo:
        int degree = 2; // Fixed here!

        JetMatrix xcj(1), ywj(1);
        Flash->molarfractions_jet(Theta2T(Theta), xcj, ywj, degree);

        xc_T         = xcj.get(0);
        dxc_dTheta   = xcj.get(0, 0)*T_typical_;
        d2xc_dTheta2 = xcj.get(0, 0, 0)*T_typical_*T_typical_;

        yw_T         = ywj.get(0);
        dyw_dTheta   = ywj.get(0, 0)*T_typical_;
        d2yw_dTheta2 = ywj.get(0, 0, 0)*T_typical_*T_typical_;

        used_Theta = Theta;
    }
    
    return;
}

//Thermodynamics::Thermodynamics(double mc, double mw, const char *hsigmaC_name){
//    liquid = new JetSinglePhaseLiquid(mc, mw, P);
//    vapor  = new JetSinglePhaseVapor(mc, mw, P);

//    info_hsigmaC = create_spline(hsigmaC_name, "hsigmaC", P, hsigmaC_);
//    printf("Thermodynamics ctor: info_hsigmaC = %d\n", info_hsigmaC);
//    if (info_hsigmaC == SPLINE_ERROR) exit(0);
//}

Thermodynamics::Thermodynamics(const char *hsigmaC_name){
    double mc = 0.044;
    double mw = 0.018;

    liquid = new JetSinglePhaseLiquid(mc, mw, P);
    vapor  = new JetSinglePhaseVapor(mc, mw, P);

    info_hsigmaC = create_spline(hsigmaC_name, "hsigmaC", P, hsigmaC_);
    printf("Thermodynamics ctor: info_hsigmaC = %d\n", info_hsigmaC);
    if (info_hsigmaC == SPLINE_ERROR) exit(0);
}

Thermodynamics::~Thermodynamics(){
    delete vapor;
    delete liquid;
}

// Generate a spline
int Thermodynamics::create_spline(const char *name, const char *verify, double P, spline1dinterpolant &spline){
    // Open the file that contains the data needed for the creation of the spline
    FILE *fid;
    fid = fopen(name, "r");
    if (fid == NULL){
        std::cout << "***** Thermodynamics::create_spline: Error while creating spline from file \"" << std::string(name) << "\" *****" << std::endl;
        return SPLINE_ERROR;
    }

    // Read the name of the variable and the pressure and verify them.
    char name_variable[100]; 
    fscanf(fid, "%s", name_variable);

    double Ptest;
    fscanf(fid, "%lf", &Ptest);
     

    if (strcmp(name_variable, verify) != 0 || Ptest != P){
        fclose(fid);

        std::cout << "***** Thermodynamics::create_spline: Error while creating spline from file \"" << std::string(name) << "\" *****" << std::endl;
        return SPLINE_ERROR;
    }

    // Read the number of base functions
    int n_base_func;
    fscanf(fid, "%d", &n_base_func);

    // Read the rest of the data
    int n;
    fscanf(fid, "%d", &n);

    ap::real_1d_array x, y;
    x.setlength(n);
    y.setlength(n);

    double xtemp, ytemp;

    for (int i = 0; i < n; i++){
        fscanf(fid, "%lf %lf", &xtemp, &ytemp);
        x(i) = xtemp;
        y(i) = ytemp;
    }

    // Close the file with the data
    fclose(fid);

    // Generate the spline
    int spline_info;
    spline1dfitreport report;
    spline1dfitcubic(x, y, n, n_base_func,
                     spline_info, 
                     spline, 
                     report);
    
    
    cout <<"Spline info: "<<spline_info<<endl;

    if (spline_info > 0) return SPLINE_OK;
    else {
        std::cout << "***** Thermodynamics::create_spline: Error while creating spline from file \"" << std::string(name) << "\" *****" << std::endl;
        return SPLINE_ERROR;
    }
}

// hsigmaC
//
int Thermodynamics::hsigmaC_jet(const double Theta, int degree, JetMatrix &hsigmaCj){
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    double T = Theta2T(Theta); // T = Theta*T_typical_ + Tref_water
    double h, d_h, d2_h;
    spline1ddiff(hsigmaC_, T, h, d_h, d2_h);

    if (degree >= 0){
        hsigmaCj.set(0, h/h_typical_);

        if (degree >= 1){
            hsigmaCj.set(0, 0, d_h*T_typical_/h_typical_);

            if (degree == 2){
                hsigmaCj.set(0, 0, 0, d2_h*T_typical_*T_typical_/h_typical_);
            }
        }
    }

    return 2; // SUCCESSFUL_PROCEDURE;
}

// Rock Enthalpy Volume
//
int Thermodynamics::RockEnthalpyVol_jet(const double Theta, int degree, JetMatrix &revj) const {
    //std::cout << "Thermodynamics::RockEnthalpyVol_jet: degree = " << degree << std::endl;

    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE

    double T = Theta2T(Theta); // T = Theta*T_typical_ + Tref_water
    //std::cout << "Thermodynamics::RockEnthalpyVol_jet: T = " << T << std::endl;

    if (degree >= 0){
        double rev = Rock_Cr*(T - Tref_rock)/(Rho_typical_*h_typical_);
//        std::cout << "Thermodynamics::RockEnthalpyVol_jet: rev = " << rev << std::endl;
//        std::cout << "Thermodynamics::RockEnthalpyVol_jet: revj.size() = " << revj.size() << std::endl;

        revj.set(0, rev);

        if (degree >= 1){
            double drev_dTheta= Rock_Cr*T_typical_/(Rho_typical_*h_typical_);

            revj.set(0, 0, drev_dTheta);

            if (degree == 2){
                double d2rev_dTheta2 = 0.0;

                revj.set(0, 0, 0, d2rev_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;
}


// Super Critical Enthalpy Volume
//
int Thermodynamics::SuperCriticEnthalpyVol_jet(const double yw, const double Theta, int degree, JetMatrix &Hsij) const {
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE

    JetMatrix rhosicj(2);
    Rhosic_jet(yw, Theta, degree, rhosicj);

    Thermodynamics *temp = const_cast<Thermodynamics*>(this);
    JetMatrix hsigmaCj(1);
    temp->hsigmaC_jet(Theta, degree, hsigmaCj);

    if (degree >= 0){
        double rhosic  = rhosicj.get(0);
        double hsigmaC = hsigmaCj.get(0);

        double Hsi = rhosic*hsigmaC;

        Hsij.set(0, Hsi);

        if (degree >= 1){
            double drhosic_dyw     = rhosicj.get(0, 0);
            double drhosic_dTheta  = rhosicj.get(0, 1);
            double dhsigmaC_dTheta = hsigmaCj.get(0, 0);

            double dHsi_dyw    = drhosic_dyw*hsigmaC;
            double dHsi_dTheta = drhosic_dTheta*hsigmaC + rhosic*dhsigmaC_dTheta;

            Hsij.set(0, 0, dHsi_dyw);
            Hsij.set(0, 1, dHsi_dTheta);

            if (degree == 2){
                double d2rhosic_dyw2      = rhosicj.get(0, 0, 0);
                double d2rhosic_dywdTheta = rhosicj.get(0, 0, 1);
                double d2rhosic_dThetadyw = rhosicj.get(0, 1, 0);
                double d2rhosic_dTheta2   = rhosicj.get(0, 1, 1);

                double d2hsigmaC_dTheta2  = hsigmaCj.get(0, 0, 0);

                double d2Hsi_dyw2      = d2rhosic_dyw2*hsigmaC;
                double d2Hsi_dywdTheta = d2rhosic_dywdTheta*hsigmaC + drhosic_dyw*dhsigmaC_dTheta;
                double d2Hsi_dThetadyw = d2Hsi_dywdTheta;
                double d2Hsi_dTheta2   = d2rhosic_dTheta2*hsigmaC + 2.0*drhosic_dTheta*dhsigmaC_dTheta +
                                         rhosic*d2hsigmaC_dTheta2;

                Hsij.set(0, 0, 0, d2Hsi_dyw2);
                Hsij.set(0, 0, 1, d2Hsi_dywdTheta);
                Hsij.set(0, 1, 0, d2Hsi_dThetadyw);
                Hsij.set(0, 1, 1, d2Hsi_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;
}

// Hsij is expected to be declared outside as:
//
//    JetMatrix Hsij(1);
//
int Thermodynamics::SuperCriticEnthalpyVol_jet(const double Theta, int degree, JetMatrix &Hsij) const {
    SetCompositions(Theta);

    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE

    JetMatrix rhosicj(1);
    Rhosic_jet(Theta, degree, rhosicj);

    Thermodynamics *temp = const_cast<Thermodynamics*>(this);
    JetMatrix hsigmaCj(1);
    temp->hsigmaC_jet(Theta, degree, hsigmaCj);

    if (degree >= 0){
        double rhosic  = rhosicj.get(0);
        double hsigmaC = hsigmaCj.get(0);

        double Hsi = rhosic*hsigmaC;

        Hsij.set(0, Hsi);

        if (degree >= 1){
            double drhosic_dTheta  = rhosicj.get(0, 0);
            double dhsigmaC_dTheta = hsigmaCj.get(0, 0);

            double dHsi_dTheta = drhosic_dTheta*hsigmaC + rhosic*dhsigmaC_dTheta;

            Hsij.set(0, 0, dHsi_dTheta);

            if (degree == 2){
                double d2rhosic_dTheta2   = rhosicj.get(0, 0, 0);

                double d2hsigmaC_dTheta2  = hsigmaCj.get(0, 0, 0);

                double d2Hsi_dTheta2   = d2rhosic_dTheta2*hsigmaC + 2.0*drhosic_dTheta*dhsigmaC_dTheta +
                                         rhosic*d2hsigmaC_dTheta2;

                Hsij.set(0, 0, 0, d2Hsi_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;

//    JetMatrix OriginalHsij(2);
//    int info = SuperCriticEnthalpyVol_jet(yw_T, Theta, degree, OriginalHsij);
//    
//    if (degree >= 0){
//        Hsij.set(0, OriginalHsij.get(0));
//        if (degree >= 1){
//            Hsij.set(0, 0, OriginalHsij.get(0, 1));
//            if (degree >= 2){
//               Hsij.set(0, 0, 0, OriginalHsij.get(0, 1, 1));
//            }
//        }
//    }

//    return info;


}
//
// End Super Critical Enthalpy Volume




// Pure water enthalpy
//
int Thermodynamics::hW_jet(const double Theta, int degree, JetMatrix &hWj){
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    double T = Theta2T(Theta); // T = Theta*T_typical_ + Tref_water

    if (degree >= 0){
        double hW = Cw_*(T - Tref_water);

        hWj.set(0, hW/h_typical_);

        if (degree >= 1){
            double dhW_dT = Cw_;

            hWj.set(0, 0, dhW_dT*T_typical_/h_typical_);

            if (degree == 2){
                double d2hW_dT2 = 0.0;

                hWj.set(0, 0, 0, d2hW_dT2*T_typical_*T_typical_/h_typical_);
            }
        }
    }

    return 2; // SUCCESSFUL_PROCEDURE;
}



// Aqueous Enthalpy Volume
//
int Thermodynamics::AqueousEnthalpyVol_jet(const double xc, const double Theta, int degree, JetMatrix &Haj) const {
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    ////printf("    AqueousEnthalpyVol: xc = %g, Theta = %g\n", xc, Theta);

    JetMatrix rhoawj(2);
    Rhoaw_jet(xc, Theta, degree, rhoawj); //printf("Done:   Rhoaw_jet\n");

//    std::cout << "AqueousEnthalpyVol_jet classic before..." << std::endl;
    Thermodynamics *temp = const_cast<Thermodynamics*>(this);
//    std::cout << "AqueousEnthalpyVol_jet classic after..." << std::endl;

    JetMatrix hWj(1);
    temp->hW_jet(Theta, degree, hWj); //printf("Done:   hW_jet\n");


    if (degree >= 0){
        double rhoaw = rhoawj.get(0);
        double hW    = hWj.get(0);

        double Ha = rhoaw*hW;

        Haj.set(0, Ha);

        if (degree >= 1){
            double drhoaw_dxc    = rhoawj.get(0, 0);
            double drhoaw_dTheta = rhoawj.get(0, 1);
            double dhW_dTheta    = hWj.get(0, 0);

            double dHa_dxc    = drhoaw_dxc*hW;
            double dHa_dTheta = drhoaw_dTheta*hW + rhoaw*dhW_dTheta;

            Haj.set(0, 0, dHa_dxc);
            Haj.set(0, 1, dHa_dTheta);

            if (degree == 2){
                double d2rhoaw_dxc2      = rhoawj.get(0, 0, 0);
                double d2rhoaw_dxcdTheta = rhoawj.get(0, 0, 1);
                double d2rhoaw_dThetadxc = rhoawj.get(0, 1, 0);
                double d2rhoaw_dTheta2   = rhoawj.get(0, 1, 1);
                double d2hW_dTheta2      = hWj.get(0, 0, 0);

                double d2Ha_dxc2      = d2rhoaw_dxc2*hW;
                double d2Ha_dxcdTheta = d2rhoaw_dxcdTheta*hW + drhoaw_dxc*dhW_dTheta;
                double d2Ha_dThetadxc = d2Ha_dxcdTheta;
                double d2Ha_dTheta2   = d2rhoaw_dTheta2*hW + 2.0*drhoaw_dTheta*dhW_dTheta + 
                                        rhoaw*d2hW_dTheta2;

                Haj.set(0, 0, 0, d2Ha_dxc2);
                Haj.set(0, 0, 1, d2Ha_dxcdTheta);
                Haj.set(0, 1, 0, d2Ha_dThetadxc);
                Haj.set(0, 1, 1, d2Ha_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;
}

// Haj is expected to be declared outside as:
//
//    JetMatrix Haj(1);
//
int Thermodynamics::AqueousEnthalpyVol_jet(const double Theta, int degree, JetMatrix &Haj) const {
//    std::cout << "AqueousEnthalpyVol_jet new before..." << std::endl;

    SetCompositions(Theta);

    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    ////printf("    AqueousEnthalpyVol: xc = %g, Theta = %g\n", xc, Theta);

    JetMatrix rhoawj(1);
    Rhoaw_jet(Theta, degree, rhoawj); //printf("Done:   Rhoaw_jet\n");

//    std::cout << "AqueousEnthalpyVol_jet classic before..." << std::endl;
    Thermodynamics *temp = const_cast<Thermodynamics*>(this);
//    std::cout << "AqueousEnthalpyVol_jet classic after..." << std::endl;

    JetMatrix hWj(1);
    temp->hW_jet(Theta, degree, hWj); //printf("Done:   hW_jet\n");

    if (degree >= 0){
        double rhoaw = rhoawj.get(0);
        double hW    = hWj.get(0);

        double Ha = rhoaw*hW;

        Haj.set(0, Ha);

        if (degree >= 1){
            double drhoaw_dTheta = rhoawj.get(0, 0);
            double dhW_dTheta    = hWj.get(0, 0);

            double dHa_dTheta = drhoaw_dTheta*hW + rhoaw*dhW_dTheta;

            Haj.set(0, 0, dHa_dTheta);

            if (degree == 2){
                double d2rhoaw_dTheta2   = rhoawj.get(0, 0, 0);
                double d2hW_dTheta2      = hWj.get(0, 0, 0);

                double d2Ha_dTheta2   = d2rhoaw_dTheta2*hW + 2.0*drhoaw_dTheta*dhW_dTheta + 
                                        rhoaw*d2hW_dTheta2;

                Haj.set(0, 0, 0, d2Ha_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;

//    JetMatrix OriginalHaj(2);
////    int info = SuperCriticEnthalpyVol_jet(xc_T, Theta, degree, OriginalHaj); 
//    int info = AqueousEnthalpyVol_jet(xc_T, Theta, degree, OriginalHaj); // Modified by Morante
////    std::cout << "AqueousEnthalpyVol_jet new after..." << std::endl;    

//    if (degree >= 0){
//        Haj.set(0, OriginalHaj.get(0));
//        if (degree >= 1){
//            Haj.set(0, 0, OriginalHaj.get(0, 1));
//            if (degree >= 2){
//               Haj.set(0, 0, 0, OriginalHaj.get(0, 1, 1));
//            }
//        }
//    }

//    return info;
}
//
// End Aqueous Enthalpy Volume

// Rho_{sigma c}
//
int Thermodynamics::Rhosic_jet(const double yw, const double Theta, int degree, JetMatrix &rhosicj) const {
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    double T = Theta2T(Theta); // T = Theta*T_typical_ + Tref_water

    JetMatrix rhosigmacj(2);
    vapor->rhosigmac_jet(yw, T, degree, rhosigmacj);

    if (degree >= 0){
        double rhosigmac = rhosigmacj.get(0);

        double rhosic = rhosigmac/Rho_typical_;
      

        rhosicj.set(0, rhosic);

        if (degree >= 1){
            double drhosigmac_dyw = rhosigmacj.get(0, 0);
            double drhosigmac_dT  = rhosigmacj.get(0, 1);

            double drhosic_dyw, drhosic_dTheta;

            drhosic_dyw = drhosigmac_dyw/Rho_typical_;

            drhosic_dTheta = drhosigmac_dT*T_typical_/Rho_typical_;

            rhosicj.set(0, 0, drhosic_dyw);
            rhosicj.set(0, 1, drhosic_dTheta);

            if (degree == 2){
                double d2rhosigmac_dyw2  = rhosigmacj.get(0, 0, 0);
                double d2rhosigmac_dywdT = rhosigmacj.get(0, 0, 1);
                double d2rhosigmac_dTdyw = rhosigmacj.get(0, 1, 0);
                double d2rhosigmac_dT2   = rhosigmacj.get(0, 1, 1);

                double d2rhosic_dyw2, d2rhosic_dywdTheta, d2rhosic_dThetadyw, d2rhosic_dTheta2;

                d2rhosic_dyw2      = d2rhosigmac_dyw2/Rho_typical_;
                d2rhosic_dywdTheta = d2rhosigmac_dywdT*T_typical_/Rho_typical_;
                d2rhosic_dThetadyw = d2rhosic_dywdTheta;
                d2rhosic_dTheta2   = d2rhosigmac_dT2*T_typical_*T_typical_/Rho_typical_;

                rhosicj.set(0, 0, 0, d2rhosic_dyw2);
                rhosicj.set(0, 0, 1, d2rhosic_dywdTheta);
                rhosicj.set(0, 1, 0, d2rhosic_dThetadyw);
                rhosicj.set(0, 1, 1, d2rhosic_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;
}

// rhosicj is expected to be declared outside as:
//
//    JetMatrix rhosicj(1);
//
int Thermodynamics::Rhosic_jet(const double Theta, int degree, JetMatrix &rhosicj) const {
    SetCompositions(Theta);

    JetMatrix Originalrhosicj(2);
    int info = Rhosic_jet(yw_T, Theta, degree, Originalrhosicj);
    
    if (degree >= 0){
        rhosicj.set(0, Originalrhosicj.get(0));
        if (degree >= 1){
            double drhosic_dyw    = Originalrhosicj.get(0, 0);
            double drhosic_dTheta = Originalrhosicj.get(0, 1);

            double dRhosic_dTheta = drhosic_dyw*dyw_dTheta + drhosic_dTheta;

            rhosicj.set(0, 0, dRhosic_dTheta);

            if (degree >= 2){
                double d2rhosic_dyw2      = Originalrhosicj.get(0, 0, 0);
                double d2rhosic_dywdTheta = Originalrhosicj.get(0, 0, 1);
                double d2rhosic_dThetadyw = Originalrhosicj.get(0, 1, 0);
                double d2rhosic_dTheta2   = Originalrhosicj.get(0, 1, 1);

                double d2Rhosic_dTheta2 = dyw_dTheta*(dyw_dTheta*d2rhosic_dyw2 + d2rhosic_dywdTheta) + 
                                          d2rhosic_dywdTheta*dyw_dTheta + drhosic_dyw*d2yw_dTheta2 + 
                                          d2rhosic_dTheta2;

                rhosicj.set(0, 0, 0, d2Rhosic_dTheta2);
            }
        }
    }

    return info;
}
//
// End Rho_{sigma c}



// Rho_{sigma w}
//
int Thermodynamics::Rhosiw_jet(const double yw, const double Theta, int degree, JetMatrix &rhosiwj) const {
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    double T = Theta2T(Theta); // T = Theta*T_typical_ + Tref_water

    JetMatrix rhosigmawj(2);
    vapor->rhosigmaw_jet(yw, T, degree, rhosigmawj);

    if (degree >= 0){
        double rhosigmaw = rhosigmawj.get(0);

        double rhosiw = rhosigmaw/Rho_typical_;

        rhosiwj.set(0, rhosiw);

        if (degree >= 1){
            double drhosigmaw_dyw = rhosigmawj.get(0, 0);
            double drhosigmaw_dT  = rhosigmawj.get(0, 1);

            double drhosiw_dyw, drhosiw_dTheta;

            drhosiw_dyw = drhosigmaw_dyw/Rho_typical_;

            drhosiw_dTheta = drhosigmaw_dT*T_typical_/Rho_typical_;
//            drhosiw_dTheta = drhosigmaw_dT;

            rhosiwj.set(0, 0, drhosiw_dyw);
            rhosiwj.set(0, 1, drhosiw_dTheta);

            if (degree == 2){
                double d2rhosigmaw_dyw2  = rhosigmawj.get(0, 0, 0);
                double d2rhosigmaw_dywdT = rhosigmawj.get(0, 0, 1);
                double d2rhosigmaw_dTdyw = rhosigmawj.get(0, 1, 0);
                double d2rhosigmaw_dT2   = rhosigmawj.get(0, 1, 1);

                double d2rhosiw_dyw2, d2rhosiw_dywdTheta, d2rhosiw_dThetadyw, d2rhosiw_dTheta2;

                d2rhosiw_dyw2      = d2rhosigmaw_dyw2/Rho_typical_;
                d2rhosiw_dywdTheta = d2rhosigmaw_dywdT*T_typical_/Rho_typical_;
                d2rhosiw_dThetadyw = d2rhosiw_dywdTheta;
                d2rhosiw_dTheta2   = d2rhosigmaw_dT2*T_typical_*T_typical_/Rho_typical_;

                rhosiwj.set(0, 0, 0, d2rhosiw_dyw2);
                rhosiwj.set(0, 0, 1, d2rhosiw_dywdTheta);
                rhosiwj.set(0, 1, 0, d2rhosiw_dThetadyw);
                rhosiwj.set(0, 1, 1, d2rhosiw_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;
}

// rhosiwj is expected to be declared outside as:
//
//    JetMatrix rhosiwj(1);
//
int Thermodynamics::Rhosiw_jet(const double Theta, int degree, JetMatrix &rhosiwj) const {
    SetCompositions(Theta);

    JetMatrix Originalrhosiwj(2);
    int info = Rhosiw_jet(yw_T, Theta, degree, Originalrhosiwj);
    
    if (degree >= 0){
        rhosiwj.set(0, Originalrhosiwj.get(0));
        if (degree >= 1){
//            rhosiwj.set(0, 0, Originalrhosiwj.get(0, 1));


            double drhosiw_dyw    = Originalrhosiwj.get(0, 0);
            double drhosiw_dTheta = Originalrhosiwj.get(0, 1);

            double dRhosiw_dTheta = drhosiw_dyw*dyw_dTheta + drhosiw_dTheta;

            rhosiwj.set(0, 0, dRhosiw_dTheta);
            if (degree >= 2){
//               rhosiwj.set(0, 0, 0, Originalrhosiwj.get(0, 1, 1));

                double d2rhosiw_dyw2      = Originalrhosiwj.get(0, 0, 0);
                double d2rhosiw_dywdTheta = Originalrhosiwj.get(0, 0, 1);
                double d2rhosiw_dThetadyw = Originalrhosiwj.get(0, 1, 0);
                double d2rhosiw_dTheta2   = Originalrhosiwj.get(0, 1, 1);

                double d2Rhosiw_dTheta2 = dyw_dTheta*(dyw_dTheta*d2rhosiw_dyw2 + d2rhosiw_dywdTheta) + 
                                          d2rhosiw_dywdTheta*dyw_dTheta + drhosiw_dyw*d2yw_dTheta2 + 
                                          d2rhosiw_dTheta2;

                rhosiwj.set(0, 0, 0, d2Rhosiw_dTheta2);
            }
        }
    }

    return info;
}
//
// End Rho_{sigma w}



// Rho_{a c}
//
int Thermodynamics::Rhoac_jet(const double xc, const double Theta, int degree, JetMatrix &rhoacj) const {
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    double T = Theta2T(Theta); // T = Theta*T_typical_ + Tref_water

    JetMatrix liquid_rhoacj(2);
    liquid->rhoac_jet(xc, T, degree, liquid_rhoacj);

    if (degree >= 0){
        double liquid_rhoac = liquid_rhoacj.get(0);

        double rhoac = liquid_rhoac/Rho_typical_;

        rhoacj.set(0, rhoac);

        if (degree >= 1){
            double dliquid_rhoac_dxc = liquid_rhoacj.get(0, 0);
            double dliquid_rhoac_dT  = liquid_rhoacj.get(0, 1);

            double drhoac_dxc, drhoac_dTheta;

            drhoac_dxc    = dliquid_rhoac_dxc/Rho_typical_;

            drhoac_dTheta = dliquid_rhoac_dT*T_typical_/Rho_typical_;

            rhoacj.set(0, 0, drhoac_dxc);
            rhoacj.set(0, 1, drhoac_dTheta);

            if (degree == 2){
                double d2liquid_rhoac_dxc2  = liquid_rhoacj.get(0, 0, 0);
                double d2liquid_rhoac_dxcdT = liquid_rhoacj.get(0, 0, 1);
                double d2liquid_rhoac_dTdxc = liquid_rhoacj.get(0, 1, 0);
                double d2liquid_rhoac_dT2   = liquid_rhoacj.get(0, 1, 1);

                double d2rhoac_dxc2, d2rhoac_dxcdTheta, d2rhoac_dThetadxc, d2rhoac_dTheta2;

                d2rhoac_dxc2      = d2liquid_rhoac_dxc2/Rho_typical_;
                d2rhoac_dxcdTheta = d2liquid_rhoac_dxcdT*T_typical_/Rho_typical_;
                d2rhoac_dThetadxc = d2rhoac_dxcdTheta;
                d2rhoac_dTheta2   = d2liquid_rhoac_dT2*T_typical_*T_typical_/Rho_typical_;

                rhoacj.set(0, 0, 0, d2rhoac_dxc2);
                rhoacj.set(0, 0, 1, d2rhoac_dxcdTheta);
                rhoacj.set(0, 1, 0, d2rhoac_dThetadxc);
                rhoacj.set(0, 1, 1, d2rhoac_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;
}

// rhoacj is expected to be declared outside as:
//
//    JetMatrix rhoacj(1);
//
int Thermodynamics::Rhoac_jet(const double Theta, int degree, JetMatrix &rhoacj) const {
    SetCompositions(Theta);

    JetMatrix Originalrhoacj(2);
    int info = Rhoac_jet(xc_T, Theta, degree, Originalrhoacj);
    
    if (degree >= 0){
        rhoacj.set(0, Originalrhoacj.get(0));
        if (degree >= 1){
//            rhoacj.set(0, 0, Originalrhoacj.get(0, 1));
            double drhoac_dxc    = Originalrhoacj.get(0, 0);
            double drhoac_dTheta = Originalrhoacj.get(0, 1);

            double dRhoac_dTheta = drhoac_dxc*dxc_dTheta + drhoac_dTheta;

            rhoacj.set(0, 0, dRhoac_dTheta);

            if (degree >= 2){
//               rhoacj.set(0, 0, 0, Originalrhoacj.get(0, 1, 1));
                double d2rhoac_dxc2      = Originalrhoacj.get(0, 0, 0);
                double d2rhoac_dxcdTheta = Originalrhoacj.get(0, 0, 1);
                double d2rhoac_dThetadxc = Originalrhoacj.get(0, 1, 0);
                double d2rhoac_dTheta2   = Originalrhoacj.get(0, 1, 1);

                double d2Rhoac_dTheta2 = dxc_dTheta*(dxc_dTheta*d2rhoac_dxc2 + d2rhoac_dxcdTheta) + 
                                         d2rhoac_dxcdTheta*dxc_dTheta + drhoac_dxc*d2xc_dTheta2 + 
                                         d2rhoac_dTheta2;

                rhoacj.set(0, 0, 0, d2Rhoac_dTheta2);

            }
        }
    }

    return info;
}
//
// End Rho_{a c}



// Rho_{a w}
//
int Thermodynamics::Rhoaw_jet(const double xc, const double Theta, int degree, JetMatrix &rhoawj) const {
    if (degree < 0 || degree > 2) return -1; // ABORTED_PROCEDURE;

    double T = Theta2T(Theta); // T = Theta*T_typical_ + Tref_water

    JetMatrix liquid_rhoawj(2); //printf("Before: liquid->rhoaw_jet\n");
    liquid->rhoaw_jet(xc, T, degree, liquid_rhoawj); //printf("Done: liquid->rhoaw_jet\n");

    if (degree >= 0){
        double liquid_rhoaw = liquid_rhoawj.get(0);

        double rhoaw = liquid_rhoaw/Rho_typical_;

        rhoawj.set(0, rhoaw);

        if (degree >= 1){
            double dliquid_rhoaw_dxc = liquid_rhoawj.get(0, 0);
            double dliquid_rhoaw_dT  = liquid_rhoawj.get(0, 1);

            double drhoaw_dxc, drhoaw_dTheta;

            drhoaw_dxc    = dliquid_rhoaw_dxc/Rho_typical_;

            drhoaw_dTheta = dliquid_rhoaw_dT*T_typical_/Rho_typical_;

            rhoawj.set(0, 0, drhoaw_dxc);
            rhoawj.set(0, 1, drhoaw_dTheta);

            if (degree == 2){
                double d2liquid_rhoaw_dxc2  = liquid_rhoawj.get(0, 0, 0);
                double d2liquid_rhoaw_dxcdT = liquid_rhoawj.get(0, 0, 1);
                double d2liquid_rhoaw_dTdxc = liquid_rhoawj.get(0, 1, 0);
                double d2liquid_rhoaw_dT2   = liquid_rhoawj.get(0, 1, 1);

                double d2rhoaw_dxc2, d2rhoaw_dxcdTheta, d2rhoaw_dThetadxc, d2rhoaw_dTheta2;

                d2rhoaw_dxc2      = d2liquid_rhoaw_dxc2/Rho_typical_;
                d2rhoaw_dxcdTheta = d2liquid_rhoaw_dxcdT*T_typical_/Rho_typical_;
                d2rhoaw_dThetadxc = d2rhoaw_dxcdTheta;
                d2rhoaw_dTheta2   = d2liquid_rhoaw_dT2*T_typical_*T_typical_/Rho_typical_;

                rhoawj.set(0, 0, 0, d2rhoaw_dxc2);
                rhoawj.set(0, 0, 1, d2rhoaw_dxcdTheta);
                rhoawj.set(0, 1, 0, d2rhoaw_dThetadxc);
                rhoawj.set(0, 1, 1, d2rhoaw_dTheta2);
            }
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE;
}

// rhoawj is expected to be declared outside as:
//
//    JetMatrix rhoawj(1);
//
int Thermodynamics::Rhoaw_jet(const double Theta, int degree, JetMatrix &rhoawj) const {
    SetCompositions(Theta);

    JetMatrix Originalrhoawj(2);
    int info = Rhoaw_jet(xc_T, Theta, degree, Originalrhoawj);
    
    if (degree >= 0){
        rhoawj.set(0, Originalrhoawj.get(0));
        if (degree >= 1){
//            rhoawj.set(0, 0, Originalrhoawj.get(0, 1));
            double drhoaw_dxc    = Originalrhoawj.get(0, 0);
            double drhoaw_dTheta = Originalrhoawj.get(0, 1);

            double dRhosic_dTheta = drhoaw_dxc*dxc_dTheta + drhoaw_dTheta;

            rhoawj.set(0, 0, dRhosic_dTheta);

            if (degree >= 2){
//               rhoawj.set(0, 0, 0, Originalrhoawj.get(0, 1, 1));
                double d2rhoaw_dxc2      = Originalrhoawj.get(0, 0, 0);
                double d2rhoaw_dxcdTheta = Originalrhoawj.get(0, 0, 1);
                double d2rhoaw_dThetadxc = Originalrhoawj.get(0, 1, 0);
                double d2rhoaw_dTheta2   = Originalrhoawj.get(0, 1, 1);

                double d2Rhoaw_dTheta2 = dxc_dTheta*(dxc_dTheta*d2rhoaw_dxc2 + d2rhoaw_dxcdTheta) + 
                                         d2rhoaw_dxcdTheta*dxc_dTheta + drhoaw_dxc*d2xc_dTheta2 + 
                                         d2rhoaw_dTheta2;

                rhoawj.set(0, 0, 0, d2Rhoaw_dTheta2);

            }
        }
    }

    return info;
}
//
// End Rho_{a w}


// Viscosity "JETS"
//
void Thermodynamics::muw(double T, double &muw, double &dmuw_dT, double &d2muw_dT2) const {

    double inv_T = 1. / T;
    double inv_T2 = inv_T*inv_T;
    double inv_T3 = inv_T2*inv_T;
    double inv_T4 = inv_T3*inv_T;
    double inv_T5 = inv_T4*inv_T;
    double inv_T6 = inv_T5*inv_T;
    double inv_T7 = inv_T6*inv_T;

    muw = (b0 + b1 * inv_T + b2 * inv_T2 + b3 * inv_T3 + b4 * inv_T4 + b5 * inv_T5);
    dmuw_dT = -(b1 * inv_T2 + 2. * b2 * inv_T3 + 3. * b3 * inv_T4 + 4. * b4 * inv_T5 + 5. * b5 * inv_T6);
    d2muw_dT2 = 2. * b1 * inv_T3 + 6. * b2 * inv_T4 + 12. * b3 * inv_T5 + 20. * b4 * inv_T6 + 30. * b5*inv_T7;

    return;
}

void Thermodynamics::inv_muw(double T, double &nuw, double &dnuw_dT, double &d2nuw_dT2) const {

    double inv_T = 1. / T;
    double inv_T2 = inv_T*inv_T;
    double inv_T3 = inv_T2*inv_T;
    double inv_T4 = inv_T3*inv_T;
    double inv_T5 = inv_T4*inv_T;
    double inv_T6 = inv_T5*inv_T;
    double inv_T7 = inv_T6*inv_T;

    double muw = (b0 + b1 * inv_T + b2 * inv_T2 + b3 * inv_T3 + b4 * inv_T4 + b5 * inv_T5);
    double dmuw_dT = -(b1 * inv_T2 + 2. * b2 * inv_T3 + 3. * b3 * inv_T4 + 4. * b4 * inv_T5 + 5. * b5 * inv_T6);
    double d2muw_dT2 = 2. * b1 * inv_T3 + 6. * b2 * inv_T4 + 12. * b3 * inv_T5 + 20. * b4 * inv_T6 + 30. * b5*inv_T7;

    nuw = 1. / muw;
    double inv_muw2 = nuw*nuw;
    double inv_muw3 = inv_muw2*nuw;

    dnuw_dT = -dmuw_dT*inv_muw2;
    d2nuw_dT2 = 2. * dmuw_dT * dmuw_dT * inv_muw3 - d2muw_dT2*inv_muw2;

    return;
}

void Thermodynamics::mug(double T, double &mug, double &dmug_dT, double &d2mug_dT2) const {

    mug = 1e-6 * (a0 * T * T * T + a1 * T * T + a2 * T + a3);
    dmug_dT = 1e-6 * (3. * a0 * T * T + 2. * a1 * T + a2);
    d2mug_dT2 = 1e-6 * (6. * a0 + 2. * a1);

    return;
}

void Thermodynamics::inv_mug(double T, double &nug, double &dnug_dT, double &d2nug_dT2) const {

    double nug_den = 1. / (a0 * T * T * T + a1 * T * T + a2 * T + a3);
    nug = 1. / (1e-6 * (a0 * T * T * T + a1 * T * T + a2 * T + a3));
    double nug_denSQ = nug_den*nug_den;

    double dmug_dT = 3. * a0 * T * T + 2. * a1 * T + a2;
    double dmug_dTSQ = dmug_dT*dmug_dT;

    dnug_dT = -1e6 * nug_denSQ*dmug_dT;
    d2nug_dT2 = 1e6 * (2. * nug_denSQ * nug_den * dmug_dTSQ - nug_denSQ * (6. * T * a0 + 2. * a1));

    return;
}


void Thermodynamics::setTtypical(double tTypical){
    T_typical_=tTypical;
}


void Thermodynamics::setRhoTypical(double rhoTypical){
    Rho_typical_=rhoTypical;
}

void Thermodynamics::UTypical(double utypical){
    U_typical_=utypical;
}

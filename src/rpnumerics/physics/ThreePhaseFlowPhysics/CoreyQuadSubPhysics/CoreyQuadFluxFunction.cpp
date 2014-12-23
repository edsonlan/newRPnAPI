#include "CoreyQuadFluxFunction.h"

CoreyQuadFluxFunction::CoreyQuadFluxFunction(Parameter *grw, Parameter *gro, Parameter *grg, 
                                             Parameter *muw, Parameter *muo, Parameter *mug,
                                             Parameter *vel){

    grw_parameter_ = grw;
    gro_parameter_ = gro;
    grg_parameter_ = grg;

    muw_parameter_ = muw;
    muo_parameter_ = muo;
    mug_parameter_ = mug;

    vel_parameter_ = vel;
}

//CoreyQuadFluxFunction::CoreyQuadFluxFunction(const CoreyQuadFluxFunction_Params &param) : FluxFunction(param) {//TODO USAR OS ELEMENTOS DO VETOR GUARDADO EM FLUXFUNCTION

//}

//CoreyQuadFluxFunction * CoreyQuadFluxFunction::clone() const {
//    return new CoreyQuadFluxFunction(*this);
//}

//CoreyQuadFluxFunction::CoreyQuadFluxFunction(const CoreyQuadFluxFunction & copy) : FluxFunction(copy.fluxParams()) {

//}

CoreyQuadFluxFunction::~CoreyQuadFluxFunction() {
}

int CoreyQuadFluxFunction::jet(const WaveState &w, JetMatrix &m, int degree) const {
    m.resize(2);

//    double grw = fluxParams().component(0);
//    double grg = fluxParams().component(1);
//    double gro = fluxParams().component(2);

//    double muw = fluxParams().component(3);
//    double mug = fluxParams().component(4);
//    double muo = fluxParams().component(5);

//    // TODO: Maybe this parameter should be always 1 (one) in the absence of gravity.
//    double vel = fluxParams().component(6);

//    double krw_p = fluxParams().component(7);
//    double krg_p = fluxParams().component(8);
//    double kro_p = fluxParams().component(9);
//
//    double cnw = fluxParams().component(10);
//    double cng = fluxParams().component(11);
//    double cno = fluxParams().component(12);


    double grw = grw_parameter_->value();
    double gro = gro_parameter_->value();
    double grg = grg_parameter_->value();

    double muw = muw_parameter_->value();
    double muo = muo_parameter_->value();
    double mug = mug_parameter_->value();

    double vel = vel_parameter_->value();

    double sw = w(0);
    double so = w(1);
    double sg = 1.0 - sw - so;
//    double CN = cnw + cng + cno;
    
    
    
    //    double CN = cnw + cng + cno;

    double kw, dkw_dsw, dkw_dso, d2kw_dsw2, d2kw_dswso, d2kw_dso2;
//    double swcnw = sw - cnw;
//    if (swcnw < 0.){ 
//        kw         = 0.;
//        dkw_dsw    = 0.;
//        dkw_dso    = 0.;

//        d2kw_dsw2  = 0.;
//        d2kw_dswso = 0.;
//        d2kw_dso2  = 0.;
//    }
//    else if (swcnw > 1. - CN) {
//        kw         = krw_p;
//        dkw_dsw    = 0.;
//        dkw_dso    = 0.;

//        d2kw_dsw2  = 0.;
//        d2kw_dswso = 0.;
//        d2kw_dso2  = 0.;
//    }
//    else {
//        double denkw = krw_p/((1. - CN)*(1. - CN));
        kw         = sw * sw; // swcnw*swcnw
        dkw_dsw    = 2. * sw;
        dkw_dso    = 0.;

        d2kw_dsw2  = 2.;
        d2kw_dswso = 0.;
        d2kw_dso2  = 0.;
//    }

    double ko, dko_dsw, dko_dso, d2ko_dsw2, d2ko_dswso, d2ko_dso2;
//    double socno = so - cno;
//    if (socno < 0.){ 
//        ko         = 0.;
//        dko_dsw    = 0.;
//        dko_dso    = 0.;

//        d2ko_dsw2  = 0.;
//        d2ko_dswso = 0.;
//        d2ko_dso2  = 0.;
//    }
//    else if (socno > 1. - CN) {
//        ko         = kro_p;
//        dko_dsw    = 0.;
//        dko_dso    = 0.;

//        d2ko_dsw2  = 0.;
//        d2ko_dswso = 0.;
//        d2ko_dso2  = 0.;
//    }
//    else {
//        double denko = kro_p/((1. - CN)*(1. - CN));
        ko         = so * so;
        dko_dsw    = 0.;
        dko_dso    = 2. * so;

        d2ko_dsw2  = 0.;
        d2ko_dswso = 0.;
        d2ko_dso2  = 2.;
//    }

    double kg, dkg_dsw, dkg_dso, d2kg_dsw2, d2kg_dswso, d2kg_dso2;
//    double sgcng = sg - cng;
//    if (sgcng < 0.){
//        kg         = 0.;
//        dkg_dsw    = 0.;
//        dkg_dso    = 0.;

//        d2kg_dsw2  = 0.;
//        d2kg_dswso = 0.;
//        d2kg_dso2  = 0.;
//    }
//    else if (sgcng > 1. - CN){
//        kg         = krg_p;
//        dkg_dsw    = 0.;
//        dkg_dso    = 0.;

//        d2kg_dsw2  = 0.;
//        d2kg_dswso = 0.;
//        d2kg_dso2  = 0.;
//    }
//    else {
//        double denkg = krg_p/((1. - CN)*(1. - CN));
        kg         = sg * sg;
        dkg_dsw    = -2. * sg;
        dkg_dso    = -2. * sg;

        d2kg_dsw2  = 2.;
        d2kg_dswso = 2.;
        d2kg_dso2  = 2.;
//    }
    
    
    
    
    
    



    

    double lkw = kw / muw; // Water mobility
    double lko = ko / muo; // Oil mobility
    double lkg = kg / mug; // Gas mobility

//        std::cout << "sg = " << sg << ", sg*sg = " << sg*sg << ", kg = " << kg << ", mug = " << mug << std::endl;

    double lk = lkw + lko + lkg; // Total mobility

    if (degree >= 0) {
        m.set(0, (lkw / lk) * (vel + lko * (grw - gro) + lkg * (grw - grg))); // fw
        m.set(1, (lko / lk) * (vel + lkw * (gro - grw) + lkg * (gro - grg))); // fo

        if (degree >= 1) {
            double ldkw_dsw = dkw_dsw / muw;
            double ldko_dsw = dko_dsw / muo;
            double ldkg_dsw = dkg_dsw / mug;
            double ldk_dsw = ldkw_dsw + ldko_dsw + ldkg_dsw;
//           std::cout << "kg_jet.get(0, 0) = " << dkg_dsw << std::endl;
//           std::cout << "ldkw_dsw = " << ldkw_dsw << ", ldko_dsw = " << ldko_dsw << ", ldkg_dsw = " << ldkg_dsw << std::endl;
//           std::cout << "ldk_dsw = " << ldk_dsw << std::endl;

            double ldkw_dso = dkw_dso / muw;
            double ldko_dso = dko_dso / muo;
            double ldkg_dso = dkg_dso / mug;
            double ldk_dso = ldkw_dso + ldko_dso + ldkg_dso;
//           std::cout << "kg_jet.get(0, 1) = " << dkg_dso << std::endl;
//           std::cout << "ldkw_dso = " << ldkw_dso << ", ldko_dso = " << ldko_dso << ", ldkg_dso = " << ldkg_dso << std::endl;
//           std::cout << "ldk_dso = " << ldk_dso << std::endl;

            double ld2kw_dsw2 = d2kw_dsw2 / muw;
            double ld2kw_dswso = d2kw_dswso / muw;
            double ld2kw_dso2 = d2kw_dso2 / muw;

            double ld2ko_dsw2 = d2ko_dsw2 / muo;
            double ld2ko_dswso = d2ko_dswso / muo;
            double ld2ko_dso2 = d2ko_dso2 / muo; // Was: d2kw_dso2 / muo;

            double ld2kg_dsw2 = d2kg_dsw2 / mug;
            double ld2kg_dswso = d2kg_dswso / mug;
            double ld2kg_dso2 = d2kg_dso2 / mug;

            // Ok.
            //
            double tw = vel + lko * (grw - gro) + lkg * (grw - grg);
            double to = vel + lkw * (gro - grw) + lkg * (gro - grg);

            // Mostly ok, but look below.
            //
            double dtodso = (gro - grg) * ldkg_dso;
            double dtodsw = (gro - grw) * ldkw_dsw + (gro - grg) * ldkg_dsw; // Was: ... + (gro - grg) * ldkg_dso <== This is an error. Marchesin & Morante, 2014-09-25, ok in Fortran.
            double dtwdso = (grw - gro) * ldko_dso + (grw - grg) * ldkg_dso;
//            double dtwdsw = ldko_dsw * (grw - gro) + (grw - grg) * ldkg_dsw; // Fortran: dtwdsw = dkgdsw * (grw - grg ), since ldko_dsw = 0.0 we use the line below. 2014-10-06.
            double dtwdsw = (grw - grg) * ldkg_dsw; // Fortran: dtwdsw = dkgdsw * (grw - grg )

            double ld2k_dsw2 = ld2kw_dsw2 + ld2ko_dsw2 + ld2kg_dsw2;
            double ld2k_dswso = ld2kw_dswso + ld2ko_dswso + ld2kg_dswso;
            double ld2k_dso2 = ld2kw_dso2 + ld2ko_dso2 + ld2kg_dso2;

            // Ok.
            //
            double zgwdsw = (lk * ldkw_dsw - lkw * ldk_dsw) / (lk * lk);
            double zgwdso = (lk * ldkw_dso - lkw * ldk_dso) / (lk * lk);
            double zgodsw = (lk * ldko_dsw - lko * ldk_dsw) / (lk * lk);
            double zgodso = (lk * ldko_dso - lko * ldk_dso) / (lk * lk);

            // Ok.
            //
            m.set(0, 0, zgwdsw * tw + (lkw / lk) * dtwdsw); // dfw_dsw
            m.set(0, 1, zgwdso * tw + (lkw / lk) * dtwdso); // dfw_dso
            m.set(1, 0, zgodsw * to + (lko / lk) * dtodsw); // dfo_dsw
            m.set(1, 1, zgodso * to + (lko / lk) * dtodso); // dfo_dso

            if (degree == 2) {
                // Ok.
                //
                double zgfww = ((lk * ld2kw_dsw2 - lkw * ld2k_dsw2) / lk - 2. * ldk_dsw * zgwdsw) / lk;
                double zgfwo = ((lk * ld2kw_dswso - lkw * ld2k_dswso) / lk - (ldk_dsw * zgwdso + zgwdsw * ldk_dso)) / lk;
                double zgfoo = ((lk * ld2kw_dso2 - lkw * ld2k_dso2) / lk - 2. * ldk_dso * zgwdso) / lk;
                double zggww = ((lk * ld2ko_dsw2 - lko * ld2k_dsw2) / lk - 2. * ldk_dsw * zgodsw) / lk;
                double zggwo = ((lk * ld2ko_dswso - lko * ld2k_dswso) / lk - (ldk_dsw * zgodso + zgodsw * ldk_dso)) / lk;
                double zggoo = ((lk * ld2ko_dso2 - lko * ld2k_dso2) / lk - 2. * ldk_dso * zgodso) / lk;

                // Ok.
                //
                double dtwdww = (grw - grg) * ld2kg_dsw2; // Was: ld2ko_dsw2 * (grw - gro) + (grw - grg) * ld2kg_dsw2;
                double dtwdwo = (grw - grg) * ld2kg_dswso; // Was: ld2ko_dswso * (grw - gro) + (grw - grg) * ld2kg_dswso;
                double dtwdoo = ld2ko_dso2 * (grw - gro) + (grw - grg) * ld2kg_dso2;
                double dtodww = ld2kw_dsw2 * (gro - grw) + (gro - grg) * ld2kg_dsw2;
                double dtodwo = (gro - grg) * ld2kg_dswso;
                double dtodoo = (gro - grg) * ld2kg_dso2;

                // Ok.
                //
                m.set(0, 0, 0, zgfww * tw + 2. * zgwdsw * dtwdsw + (lkw / lk) * dtwdww); // d2fw_dsw2;
                m.set(0, 0, 1, zgfwo * tw + zgwdsw * dtwdso + zgwdso * dtwdsw + (lkw / lk) * dtwdwo); // d2fw_dswso;
                m.set(0, 1, 0, m.get(0, 0, 1)); // d2fw_dsosw;
                m.set(0, 1, 1, zgfoo * tw + 2. * zgwdso * dtwdso + (lkw / lk) * dtwdoo); // d2fw_dso2;

                m.set(1, 0, 0, zggww * to + 2. * zgodsw * dtodsw + (lko / lk) * dtodww); // d2fo_dsw2;
                m.set(1, 0, 1, zggwo * to + zgodso * dtodsw + zgodsw * dtodso + (lko / lk) * dtodwo); // d2fo_dswso;
                m.set(1, 1, 0, m.get(1, 0, 1)); // d2fo_dsosw;
                m.set(1, 1, 1, zggoo * to + 2. * zgodso * dtodso + (lko / lk) * dtodoo); // d2fo_dso2;
            }
        }
    }
    return 2; //SUCCESSFUL_PROCEDURE;
}


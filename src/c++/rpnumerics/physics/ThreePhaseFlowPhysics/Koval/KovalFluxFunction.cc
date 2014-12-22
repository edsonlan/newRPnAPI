#include "KovalFluxFunction.h"

KovalFluxFunction::KovalFluxFunction(Parameter *grw, Parameter *gro, Parameter *grg, 
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

KovalFluxFunction::~KovalFluxFunction(){
}

int KovalFluxFunction::jet(const WaveState &w, JetMatrix &m, int degree) const {
    m.resize(2);

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

    if (degree >= 0) {
        double kw = sw;
        double ko = so;
        double kg = 1.0 - sw - so;

        double lkw = kw / muw; // Water mobility
        double lko = ko / muo; // Oil mobility
        double lkg = kg / mug; // Gas mobility

        double lk = lkw + lko + lkg; // Total mobility

        m.set(0, (lkw / lk) * (vel + lko * (grw - gro) + lkg * (grw - grg))); // fw
        m.set(1, (lko / lk) * (vel + lkw * (gro - grw) + lkg * (gro - grg))); // fo

        if (degree >= 1) {
            double dkw_dsw = 1.0;
            double dko_dsw = 0.0;
            double dkg_dsw = -1.0;

            double ldkw_dsw = dkw_dsw / muw;
            double ldko_dsw = dko_dsw / muo;
            double ldkg_dsw = dkg_dsw / mug;
            double ldk_dsw = ldkw_dsw + ldko_dsw + ldkg_dsw;

            double dkw_dso = 0.0;
            double dko_dso = 1.0;
            double dkg_dso = -1.0;

            double ldkw_dso = dkw_dso / muw;
            double ldko_dso = dko_dso / muo;
            double ldkg_dso = dkg_dso / mug;
            double ldk_dso = ldkw_dso + ldko_dso + ldkg_dso;

            double to = vel + lkw * (gro - grw) + lkg * (gro - grg);
            double tw = vel + lko * (grw - gro) + lkg * (grw - grg);

            double dtodso = (gro - grg) * ldkg_dso;
            double dtodsw = (gro - grw) * ldkw_dsw + (gro - grg) * ldkg_dsw;
            double dtwdso = (grw - gro) * ldko_dso + (grw - grg) * ldkg_dso;
            double dtwdsw = (grw - grg) * ldkg_dsw;

            double zgwdsw = (lk * ldkw_dsw - lkw * ldk_dsw) / (lk * lk);
            double zgwdso = (lk * ldkw_dso - lkw * ldk_dso) / (lk * lk);
            double zgodsw = (lk * ldko_dsw - lko * ldk_dsw) / (lk * lk);
            double zgodso = (lk * ldko_dso - lko * ldk_dso) / (lk * lk);

            m.set(0, 0, zgwdsw * tw + (lkw / lk) * dtwdsw); // dfw_dsw
            m.set(0, 1, zgwdso * tw + (lkw / lk) * dtwdso); // dfw_dso
            m.set(1, 0, zgodsw * to + (lko / lk) * dtodsw); // dfo_dsw
            m.set(1, 1, zgodso * to + (lko / lk) * dtodso); // dfo_dso

            if (degree == 2) {
//                double d2kw_dsw2  = 0.0;
//                double d2kw_dswso = 0.0;
//                double d2kw_dso2  = 0.0;

//                double ld2kw_dsw2  = 0.0;
//                double ld2kw_dswso = 0.0;
//                double ld2kw_dsosw = 0.0;
//                double ld2kw_dso2  = 0.0;

//                double ld2ko_dsw2  = 0.0;
//                double ld2ko_dswso = 0.0;
//                double ld2ko_dsosw = 0.0;
//                double ld2ko_dso2  = 0.0;

//                double ld2kg_dsw2  = 0.0;
//                double ld2kg_dswso = 0.0;
//                double ld2kg_dsosw = 0.0;
//                double ld2kg_dso2  = 0.0;

//                double ld2k_dsw2  = 0.0;
//                double ld2k_dswso = 0.0;
//                double ld2k_dso2  = 0.0;

                double zgfww = -2.*ldk_dsw*zgwdsw/lk;
                double zgfwo = -(ldk_dsw * zgwdso + zgwdsw * ldk_dso) / lk;
                double zgfoo = -2.*ldk_dso*zgwdso/lk;
                double zggww = -2.*ldk_dsw*zgodsw/lk;
                double zggwo =  -(ldk_dsw*zgodso + zgodsw*ldk_dso)/lk;
                double zggoo = -2.*ldk_dso*zgodso/lk;

//                // These terms are all NULL
//                double dtwdww = ld2ko_dsw2 * (grw - gro) + (grw - grg) * ld2kg_dsw2;
//                double dtwdwo = ld2ko_dswso * (grw - gro) + (grw - grg) * ld2kg_dswso;
//                double dtwdoo = ld2ko_dso2 * (grw - gro) + (grw - grg) * ld2kg_dso2;
//                double dtodww = ld2kw_dsw2 * (gro - grw) + (gro - grg) * ld2kg_dsw2;
//                double dtodwo = (gro - grg) * ld2kg_dswso;
//                double dtodoo = (gro - grg) * ld2kg_dso2;

                m.set(0, 0, 0, zgfww * tw + 2. * zgwdsw * dtwdsw); // d2fw_dsw2;
                m.set(0, 0, 1, zgfwo * tw + zgwdsw * dtwdso + zgwdso * dtwdsw); // d2fw_dswso;
                m.set(0, 1, 0, zgfwo * tw + zgwdsw * dtwdso + zgwdso * dtwdsw); // d2fw_dsosw;
                m.set(0, 1, 1, zgfoo * tw + 2. * zgwdso * dtwdso); // d2fw_dso2;

                m.set(1, 0, 0, zggww * to + 2. * zgodsw * dtodsw); // d2fo_dsw2;
                m.set(1, 0, 1, zggwo * to + zgodso * dtodsw + zgodsw * dtodso); // d2fo_dswso;
                m.set(1, 1, 0, zggwo * to + zgodso * dtodsw + zgodsw * dtodso); // d2fo_dsosw;
                m.set(1, 1, 1, zggoo * to + 2. * zgodso * dtodso); // d2fo_dso2;
            }
        }
    }
    return 2; //SUCCESSFUL_PROCEDURE;
}


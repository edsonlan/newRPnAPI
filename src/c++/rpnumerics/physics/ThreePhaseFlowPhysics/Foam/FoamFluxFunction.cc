#include "FoamFluxFunction.h"

FoamFluxFunction::FoamFluxFunction(Parameter *grw, Parameter *grg, Parameter *gro,
                                   Parameter *muw, Parameter *muo,
                                   Parameter *vel,
                                   FoamPermeability *fp,
                                   FoamViscosity *fv) : ThreePhaseFlowFluxFunction(){
    grw_ = grw;
    grg_ = grg;
    gro_ = gro;

    muw_ = muw;
    muo_ = muo;

    vel_ = vel;

    permeability_     = fp;
    foamvisc_ = fv;
}

FoamFluxFunction::~FoamFluxFunction() {
}

int FoamFluxFunction::jet(const RealVector &w, JetMatrix &m, int degree) const {
    m.resize(2);

    double grw = grw_->value();
    double gro = gro_->value();
    double grg = grg_->value();

    double muw = muw_->value();
    double muo = muo_->value();

    double vel = vel_->value();

    double sw = w(0);
    double so = w(1);
    double sg = 1.0 - sw - so;

    JetMatrix kw_jet;
    permeability_->PermeabilityWater_jet(w, degree, kw_jet);

    JetMatrix ko_jet;
    permeability_->PermeabilityOil_jet(w, degree, ko_jet);

    JetMatrix kg_jet;
    permeability_->PermeabilityGas_jet(w, degree, kg_jet);

    JetMatrix mug_jet;
    foamvisc_->gas_viscosity_jet(w, degree, mug_jet);

    if (degree >= 0){
        // Water mobility.
        //
        double kw = kw_jet.get(0);
        double lkw = kw/muw;

        // Oil mobility.
        //
        double ko = ko_jet.get(0);
        double lko = ko/muo;

        // Gas mobility.
        //
        double mug = mug_jet.get(0);
        double kg = kg_jet.get(0);
        double lkg = kg/mug; // ***

        // Total mobility.
        //
        double lk = lkw + lko + lkg;

        m.set(0, ( lkw/lk ) * ( vel + lko*(grw - gro) + lkg*(grw - grg) )); // fw
        m.set(1, ( lko/lk ) * ( vel + lkw*(gro - grw) + lkg*(gro - grg) )); // fo

        if (degree >= 1){ 
           double ldkw_dsw = kw_jet.get(0, 0)/muw; // dkw_dsw = kw_jet.get(0, 0)
           double ldko_dsw = ko_jet.get(0, 0)/muo; // dko_dsw = ko_jet.get(0, 0)
           double ldkg_dsw = kg_jet.get(0, 0)/mug  -  lkg*mug_jet.get(0, 0)/mug; // dkg_dsw = kg_jet.get(0, 0), dmug_dsw = mug_jet.get(0, 0) ***
           double ldk_dsw  = ldkw_dsw + ldko_dsw + ldkg_dsw;

           double ldkw_dso = kw_jet.get(0, 1)/muw; // dkw_dsw = kw_jet.get(0, 0)
           double ldko_dso = ko_jet.get(0, 1)/muo; // dko_dsw = ko_jet.get(0, 0)
           double ldkg_dso = kg_jet.get(0, 1)/mug  -  lkg*mug_jet.get(0, 1)/mug; // dkg_dsw = kg_jet.get(0, 0), dmug_dsw = mug_jet.get(0, 0) ***
           double ldk_dso  = ldkw_dso + ldko_dso + ldkg_dso;

           double tw = vel + lko*(grw - gro) + lkg*(grw - grg);
           double to = vel + lkw*(gro - grw) + lkg*(gro - grg);

           double dtodso = (gro - grg)*ldkg_dso;
           double dtodsw = (gro - grw)*ldkw_dsw + (gro - grg)*ldkg_dsw;
           double dtwdso = (grw - gro)*ldko_dso + (grw - grg)*ldkg_dso;
           double dtwdsw = (grw - gro)*ldko_dsw;// + (grw - grg)*ldkg_dsw;

           double zgwdsw = ( lk * ldkw_dsw  -  lkw * ldk_dsw ) / (lk * lk);
           double zgwdso = ( lk * ldkw_dso  -  lkw * ldk_dso ) / (lk * lk);
           double zgodsw = ( lk * ldko_dsw  -  lko * ldk_dsw ) / (lk * lk);
           double zgodso = ( lk * ldko_dso  -  lko * ldk_dso ) / (lk * lk);

           m.set(0, 0, zgwdsw*tw + (lkw/lk)*dtwdsw); // dfw_dsw
           m.set(0, 1, zgwdso*tw + (lkw/lk)*dtwdso); // dfw_dso
           m.set(1, 0, zgodsw*to + (lko/lk)*dtodsw); // dfo_dsw
           m.set(1, 1, zgodso*to + (lko/lk)*dtodso); // dfo_dso

           if (degree == 2){
               double ld2kw_dsw2  = kw_jet.get(0, 0, 0)/muw; // d2kw_dsw2  = kw_jet.get(0, 0, 0)
               double ld2kw_dswso = kw_jet.get(0, 0, 1)/muw; // d2kw_dswso = kw_jet.get(0, 0, 1)
               double ld2kw_dso2  = kw_jet.get(0, 1, 1)/muw; // d2kw_dso2  = kw_jet.get(0, 1, 1)

               double ld2ko_dsw2  = ko_jet.get(0, 0, 0)/muo; // d2kw_dsw2  = kw_jet.get(0, 0, 0)
               double ld2ko_dswso = ko_jet.get(0, 0, 1)/muo; // d2kw_dswso = kw_jet.get(0, 0, 1)
               double ld2ko_dso2  = ko_jet.get(0, 1, 1)/muo; // d2kw_dso2  = kw_jet.get(0, 1, 1)

               double ld2kg_dsw2  = kg_jet.get(0, 0, 0)/mug - kg_jet.get(0, 0)*mug_jet.get(0, 0)/(mug*mug) - 
                                    (ldkg_dsw*mug_jet.get(0, 0) + lkg*mug_jet.get(0, 0, 0))/mug + 
                                    (lkg*mug_jet.get(0, 0)*mug_jet.get(0, 0))/(mug*mug); // ***

               double ld2kg_dswso = kg_jet.get(0, 0, 1)/mug - kg_jet.get(0, 0)*mug_jet.get(0, 1)/(mug*mug) - 
                                    (ldkg_dso*mug_jet.get(0, 0) + lkg*mug_jet.get(0, 0, 1))/mug + 
                                    (lkg*mug_jet.get(0, 0)*mug_jet.get(0, 1))/(mug*mug); // ***

               double ld2kg_dso2  = kg_jet.get(0, 1, 1)/mug - kg_jet.get(0, 1)*mug_jet.get(0, 1)/(mug*mug) - 
                                    (ldkg_dso*mug_jet.get(0, 1) + lkg*mug_jet.get(0, 1, 1))/mug + 
                                    (lkg*mug_jet.get(0, 1)*mug_jet.get(0, 1))/(mug*mug); // ***

               double ld2k_dsw2  = ld2kw_dsw2 + ld2ko_dsw2 + ld2kg_dsw2; 
               double ld2k_dswso = ld2kw_dswso + ld2ko_dswso + ld2kg_dswso;
               double ld2k_dso2  = ld2kw_dso2 + ld2ko_dso2 + ld2kg_dso2;

               double zgfww =( ( lk * ld2kw_dsw2 - lkw * ld2k_dsw2 )   / lk  -  2. * ldk_dsw * zgwdsw )/ lk;
               double zgfwo =( ( lk * ld2kw_dswso - lkw * ld2k_dswso ) / lk  -  (ldk_dsw * zgwdso + zgwdsw * ldk_dso) ) / lk;
               double zgfoo =( ( lk * ld2kw_dso2 - lkw * ld2k_dso2 )   / lk  -  2. * ldk_dso * zgwdso )/ lk;
               double zggww =( ( lk * ld2ko_dsw2 - lko * ld2k_dsw2 )   / lk  -  2. * ldk_dsw * zgodsw )/ lk;
               double zggwo =( ( lk * ld2ko_dswso - lko * ld2k_dswso ) / lk  -  (ldk_dsw * zgodso + zgodsw * ldk_dso) ) / lk;
               double zggoo =( ( lk * ld2ko_dso2 - lko * ld2k_dso2 )   / lk  -  2. * ldk_dso * zgodso )/ lk;
          
               double dtwdww = ld2ko_dsw2*(grw - gro) + (grw - grg)*ld2kg_dsw2;
               double dtwdwo = ld2ko_dswso*(grw - gro) + (grw - grg)*ld2kg_dswso;
               double dtwdoo = ld2ko_dso2*(grw - gro) + (grw - grg)*ld2kg_dso2;
               double dtodww = ld2kw_dsw2*(gro - grw) + (gro - grg)*ld2kg_dsw2;
               double dtodwo = (gro - grg)*ld2kg_dswso;
               double dtodoo = (gro - grg)*ld2kg_dso2;

               m.set(0, 0, 0, zgfww*tw + 2.*zgwdsw*dtwdsw + (lkw/lk)*dtwdww); // d2fw_dsw2;
               m.set(0, 0, 1, zgfwo*tw + zgwdsw*dtwdso + zgwdso*dtwdsw + (lkw/lk)*dtwdwo); // d2fw_dswso;
               m.set(0, 1, 0, m.get(0, 0, 1));                                // d2fw_dsosw;
               m.set(0, 1, 1, zgfoo*tw + 2.*zgwdso*dtwdso + (lkw/lk)*dtwdoo); // d2fw_dso2;

               m.set(1, 0, 0, zggww*to + 2.*zgodsw*dtodsw + (lko/lk)*dtodww); // d2fo_dsw2;
               m.set(1, 0, 1, zggwo*to + zgodso*dtodsw + zgodsw*dtodso + (lko/lk)*dtodwo); // d2fo_dswso;
               m.set(1, 1, 0, m.get(1, 0, 1)); // d2fo_dsosw;
               m.set(1, 1, 1, zggoo*to + 2.*zgodso*dtodso + (lko/lk)*dtodoo); // d2fo_dso2;

            }
        }
        else return 0;
    }
    return 2; //SUCCESSFUL_PROCEDURE;
}


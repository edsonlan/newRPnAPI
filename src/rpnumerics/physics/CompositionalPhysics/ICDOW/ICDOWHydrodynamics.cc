#include "ICDOWHydrodynamics.h"

ICDOWHydrodynamics::ICDOWHydrodynamics(Parameter *muw, Parameter *muo,
                                       Parameter *grw, Parameter *gro,
                                       Parameter *vel,
                                       Parameter *swc, Parameter *lambda)
                                       
{
    muw_parameter    = muw;
    muo_parameter    = muo;
    grw_parameter    = grw;
    gro_parameter    = gro;
    vel_parameter    = vel;
    swc_parameter    = swc;
    lambda_parameter = lambda;
}

ICDOWHydrodynamics::~ICDOWHydrodynamics(){
}

void ICDOWHydrodynamics::water_fractional_flow(double sw, int degree, JetMatrix &fw_jet){
    double muw    = muw_parameter->value();
    double muo    = muo_parameter->value();
    double grw    = grw_parameter->value();
    double gro    = gro_parameter->value();
    double vel    = vel_parameter->value();
    double swc    = swc_parameter->value();
    double lambda = lambda_parameter->value();

    fw_jet.resize(1);

    if (degree >= 0){
        double so = 1.0 - sw;
 
        double kw = sw*sw;
        double ko = so*so;

        double lkw = kw / muw; // Water mobility
        double lko = ko / muo; // Oil mobility
        double lk = lkw + lko;
        double inv_lk = 1.0/lk;

        double t = lkw*inv_lk;

        double fw = t*(vel + lko*(grw - gro));
        fw_jet.set(0, fw);

        if (degree >= 1){
            double dkw_dsw  = 2.0*sw;
            double dko_dsw  = -2.0*so;
            double dlkw_dsw = dkw_dsw/muw;
            double dlko_dsw = dko_dsw/muo;
            double dlk_dsw  = dlkw_dsw + dlko_dsw;
            double dt_dsw   = (dlkw_dsw - t*dlk_dsw)*inv_lk;

            double dfw_dsw = dt_dsw*(vel + lko*(grw - gro)) + t*dlko_dsw*(grw - gro);
            fw_jet.set(0, 0, dfw_dsw);

            if (degree >= 2){
                double d2kw_dsw2  = 2.0;
                double d2ko_dsw2  = 2.0;
                double d2lkw_dsw2 = d2kw_dsw2/muw;
                double d2lko_dsw2 = d2ko_dsw2/muo;
                double d2lk_dsw2  = d2lkw_dsw2 + d2lko_dsw2;
                double d2t_dsw2   = inv_lk*(d2lkw_dsw2 - (dt_dsw*dlk_dsw + t*d2lk_dsw2) + inv_lk*(t*dlk_dsw*dlk_dsw - dlkw_dsw*dlk_dsw));

                double d2fw_dsw2 = d2t_dsw2*(vel + lko*(grw - gro)) + dt_dsw*dlko_dsw*(grw - gro) +
                                   dt_dsw*dlko_dsw*(grw - gro) + t*d2lko_dsw2*(grw - gro);
                fw_jet.set(0, 0, 0, d2fw_dsw2);
            }
        }
    }

    return;
}

//void ICDOWHydrodynamics::water_fractional_flow(double sw, int degree, JetMatrix &fw_jet){
//    double muw    = muw_parameter->value();
//    double muo    = muo_parameter->value();
//    double grw    = grw_parameter->value();
//    double gro    = gro_parameter->value();
//    double vel    = vel_parameter->value();
//    double swc    = swc_parameter->value();
//    double lambda = lambda_parameter->value();

//    fw_jet.resize(1);

//    if (degree >= 0){
//        double se;
//        if (sw > swc) se = (sw - swc)/(1.0 - swc);
//        else          se = 0.0;

//        double se2         = se*se;
//        double oneminusse  = (1.0 - se);
//        double oneminusse2 = oneminusse*oneminusse;

//        double KE = pow(se, 2.0/lambda);
//        double krw = KE*pow(se, 3.0);
//        double kro= (1.0 - KE*se)*oneminusse2;
//        double lamw=krw/muw;
//        double lamo=kro/muo;
//        double Fw_n=lamw/(lamw + lamo);
//        double fw=Fw_n*(vel + kro*(grw - gro));

//        fw_jet.set(0, fw);

//        if (degree >= 1){
//            double Dlamw= -((2.0/lambda+3.0)*se*se*KE)/(muw*(swc-1.0));
//            double Dlamo=((2/lambda+1)*oneminusse2*KE/(muo*(swc-1.0)))-(2.0*(1.0-KE*se)*oneminusse/(muo*(swc-1.0)));
//            double Dkrw=-(2/lambda+3)*KE/(se2*(swc-1));
//            double Dkro=((2/lambda+1)*oneminusse2*KE/(swc-1.0))-(2.0*(1.0-KE*se)*oneminusse/(swc-1.0));
//            double Dfw1=(Dlamw*lamo - Dlamo*lamw)/pow(lamw + lamo, 2.0);
//            double Dfw=Dfw1*(vel+kro*(grw-gro))+Fw_n*Dkro*(grw-gro);
//            fw_jet.set(0, 0, Dfw);

//            if (degree >= 2){
//                double DDlamw= ((2.0/lambda+2.0)*(2.0/lambda+3.0)*se*KE)/(muw*pow((swc-1.0),2.0));
//                double DDlamo= ((4.0*(2.0/lambda+1.0)*(1.0-se)*KE)/muo*pow(swc-1.0,2.0))-
//                               (2.0*(((se*KE)-1.0)/muo*pow(swc-1.0,2.0)))-
//                               ((2.0/lambda+1)*pow((1.0-se),2.0)*(KE)/lambda*muo*se*pow((swc-1.0),2.0));

//                double DDkro= (4*(2/lambda+1)*(1-se)*KE)-(2*(1-KE*se))/(pow((swc-1),2.0))-
//                              (2*(2/lambda+1)*(pow((1-se),2.0)*KE)/(se*lambda*pow((swc-1),2.0)));

//                double DDkrw= (2/lambda+2)*(2/lambda+3)*KE*se/pow((swc-1),2.0);

//                double DDfw1 = ((lamo*DDlamw-lamw*DDlamo)/pow(lamw+lamo, 2.0))-((2*Dfw1*(Dlamw+Dlamo))/(lamw+lamo));
//                double DDfw=(DDfw1*(vel+kro*(grw-gro))+Dfw1*Dkro*(grw-gro))+(Dfw1*Dkro+Fw_n*DDkro)*(grw-gro);

//                fw_jet.set(0, 0, 0, DDfw);
//            }
//        }
//    }

//    return;
//}

void ICDOWHydrodynamics::oil_fractional_flow(double sw, int degree, JetMatrix &fo_jet){
    double muw    = muw_parameter->value();
    double muo    = muo_parameter->value();
    double grw    = grw_parameter->value();
    double gro    = gro_parameter->value();
    double vel    = vel_parameter->value();
    double swc    = swc_parameter->value();
    double lambda = lambda_parameter->value();

    fo_jet.resize(1);

    if (degree >= 0){
        double so = 1.0 - sw;
 
        double kw = sw*sw;
        double ko = so*so;

        double lkw = kw / muw; // Water mobility
        double lko = ko / muo; // Oil mobility
        double lk = lkw + lko;
        double inv_lk = 1.0/lk;

        double to = lko*inv_lk;

        double fo = to*(vel + lkw*(gro - grw));
        fo_jet.set(0, fo);

        if (degree >= 1){
            double dkw_dsw  = 2.0*sw;
            double dko_dsw  = -2.0*so;
            double dlkw_dsw = dkw_dsw/muw;
            double dlko_dsw = dko_dsw/muo;
            double dlk_dsw  = dlkw_dsw + dlko_dsw;
            double dto_dsw   = (dlko_dsw - to*dlk_dsw)*inv_lk;

            double dfo_dsw = dto_dsw*(vel + lkw*(gro - grw)) + to*dlkw_dsw*(gro - grw);

            fo_jet.set(0, 0, dfo_dsw);

            if (degree >= 2){
                double d2kw_dsw2  = 2.0;
                double d2ko_dsw2  = 2.0;
                double d2lkw_dsw2 = d2kw_dsw2/muw;
                double d2lko_dsw2 = d2ko_dsw2/muo;
                double d2lk_dsw2  = d2lkw_dsw2 + d2lko_dsw2;
                double d2to_dsw2   = inv_lk*(d2lko_dsw2 - (dto_dsw*dlk_dsw + to*d2lk_dsw2) + inv_lk*(to*dlk_dsw*dlk_dsw - dlko_dsw*dlk_dsw));

                double d2fo_dsw2 = d2to_dsw2*(vel + lkw*(gro - grw)) + dto_dsw*dlkw_dsw*(gro - grw) +
                                   dto_dsw*dlkw_dsw*(gro - grw) + to*d2lkw_dsw2*(gro - grw);
                fo_jet.set(0, 0, 0, d2fo_dsw2);
            }
        }
    }

    return;
}

//void ICDOWHydrodynamics::oil_fractional_flow(double sw, int degree, JetMatrix &fo_jet){
//    double muw    = muw_parameter->value();
//    double muo    = muo_parameter->value();
//    double grw    = grw_parameter->value();
//    double gro    = gro_parameter->value();
//    double vel    = vel_parameter->value();
//    double swc    = swc_parameter->value();
//    double lambda = lambda_parameter->value();

//    fo_jet.resize(1);

//    if (degree >= 0){
//        double se;
//        if (sw > swc) se = (sw - swc)/(1.0 - swc);
//        else          se = 0.0;

//        double KE = pow(-(sw-swc)/(swc-1.0),2.0/lambda);
//        double krw=pow(se,3.0)*KE;
//        double kro= (1.0-KE*se)*pow(1.0 - se,2.0);
//        double lamw=krw/muw;
//        double lamo=kro/muo;
//        double Fw_n=(lamw/(lamw+lamo));
//        double fo=Fw_n*(vel+krw*(gro-grw));

//        fo_jet.set(0, fo);

//        if (degree >= 1){
//            double Dlamw= -((2.0/lambda+3.0)*pow(se,2.0)*KE)/(muw*(swc-1.0));
//            double Dlamo=(((1.0-se)*(2.0/lambda+1.0)*KE)/(muo*(swc-1.0)))-((KE*se-1.0)/(muo*(1.0-swc)));
//            double Dkrw=-(2/lambda+3)*KE/(pow(se,2.0)*(swc-1));
//            double Dkro=((2/lambda+1)*pow((1.0-se),2.0)*KE/(swc-1.0))-(2.0*(1.0-KE*se)*(1.0-se)/(swc-1.0));
//            double Dfw1=(Dlamw*(lamo)-Dlamo*(lamw))/pow((lamw+lamo),2.0);
//            double Dfo=Dfw1*(vel+krw*(gro-grw))+Fw_n*Dkrw*(gro-grw);

//            fo_jet.set(0, 0, Dfo);

//            if (degree >= 2){
//                double DDlamw= ((2.0/lambda+2.0)*(2.0/lambda+3.0)*se*KE)/(muw*pow((swc-1.0),2.0));
//                double DDlamo=((4.0*(2.0/lambda+1.0)*(1.0-se)*KE)/muo*(swc-1.0))-
//(2.0*(((se*KE)-1.0)/muo*pow(swc-1.0,2.0)))-
//((2.0/lambda+1)*pow((1.0-se),2.0)*(KE)/lambda*muo*se*pow((swc-1.0),2.0));
//                double DDkro= (4*(2/lambda+1)*(1-se)*KE)-(2*(1-KE*se))/(pow((swc-1),2.0))-
//      (2*(2/lambda+1)*(pow((1-se),2.0)*KE)/(se*lambda*pow((swc-1),2.0)));
//                double DDkrw= (2/lambda+2)*(2/lambda+3)*KE*se/pow((swc-1),2.0);

//                double DDfw1=((lamo*DDlamw-lamw*DDlamo)/pow((lamw+lamo),2.0))-2*Dfw1*(Dlamw+Dlamo)*(lamw+lamo);
//                double DDfo=(DDfw1*(vel+krw*(gro-grw)))+Dfw1*Dkrw*(grw-gro)+(Dfw1*Dkrw+Fw_n*DDkrw)*(gro-grw);

//                fo_jet.set(0, 0, 0, DDfo);
//            }
//        }
//    }

//    return;
//}



#include "Brooks_CoreyFluxFunction.h"

Brooks_CoreyFluxFunction::Brooks_CoreyFluxFunction(Parameter *muw, Parameter *muo, Parameter *mug, Brooks_CoreyPermeability *p) : FluxFunction(){
    muw_parameter = muw;
    muo_parameter = muo;
    mug_parameter = mug;

    perm = p;
}

Brooks_CoreyFluxFunction::~Brooks_CoreyFluxFunction(){
}

int Brooks_CoreyFluxFunction::jet(const WaveState &u, JetMatrix &f, int degree) const {
    // std::cout << "Jet begins." << std::endl;

    f.resize(2);

    RealVector state(3);
    state(0) = u(0);
    state(1) = u(1);
    state(2) = 1.0 - u(0) - u(1);

    // std::cout << "state = " << state << std::endl;

    JetMatrix water;
    perm->PermeabilityWater_jet(state, degree, water);

    JetMatrix gas;
    perm->PermeabilityGas_jet(state, degree, gas);

    JetMatrix oil;
    perm->PermeabilityOil_jet(state, degree, oil);

    if (degree >= 0){
        double muw = muw_parameter->value();
        double muo = muo_parameter->value();
        double mug = mug_parameter->value();

        double kw = water.get(0);
        double ko = oil.get(0);
        double kg = gas.get(0);

        double lambdaw = kw/muw;
        double lambdao = ko/muo;
        double lambdag = kg/mug;

        double lambda = lambdaw + lambdao + lambdag;
        double inv_lambda = 1.0/lambda;

        double fw = (lambdaw/lambda);//*(vel + lambdao * (grw - gro) + lambdag * (grw - grg)));
        double fo = (lambdao/lambda);//*(vel + lambdaw * (gro - grw) + lambdag * (gro - grg)));

        f.set(0, fw);
        f.set(1, fo);

        if (degree >= 1){
            double inv_lambda2 = inv_lambda*inv_lambda;

            double dlambdaw_dsw = water.get(0, 0)/muw;
            double dlambdao_dsw = oil.get(0, 0)/muo;
            double dlambdag_dsw = gas.get(0, 0)/mug;
            double dlambda_dsw  = dlambdaw_dsw + dlambdao_dsw + dlambdag_dsw;

            double dlambdaw_dso = water.get(0, 1)/muw;
            double dlambdao_dso = oil.get(0, 1)/muo;
            double dlambdag_dso = gas.get(0, 1)/mug;
            double dlambda_dso  = dlambdaw_dso + dlambdao_dso + dlambdag_dso;

            double dfw_dsw = dlambdaw_dsw/lambda - dlambda_dsw*lambdaw*inv_lambda2;
            double dfw_dso = dlambdaw_dso/lambda - dlambda_dso*lambdaw*inv_lambda2;
            double dfo_dsw = dlambdao_dsw/lambda - dlambda_dsw*lambdao*inv_lambda2;
            double dfo_dso = dlambdao_dso/lambda - dlambda_dso*lambdao*inv_lambda2;

            f.set(0, 0, dfw_dsw);
            f.set(0, 1, dfw_dso);
            f.set(1, 0, dfo_dsw);
            f.set(1, 1, dfo_dso);

            if (degree >= 2){
                double twice_inv_lambda3 = 2.0*inv_lambda2*inv_lambda;

                double d2lambdaw_dsw2  = water.get(0, 0, 0);
                double d2lambdaw_dswso = water.get(0, 0, 1);
                double d2lambdaw_dsosw = water.get(0, 1, 0);
                double d2lambdaw_dso2  = water.get(0, 1, 1);

                double d2lambdao_dsw2  = oil.get(0, 0, 0);
                double d2lambdao_dswso = oil.get(0, 0, 1);
                double d2lambdao_dsosw = oil.get(0, 1, 0);
                double d2lambdao_dso2  = oil.get(0, 1, 1);

                double d2lambdag_dsw2  = gas.get(0, 0, 0);
                double d2lambdag_dswso = gas.get(0, 0, 1);
                double d2lambdag_dsosw = gas.get(0, 1, 0);
                double d2lambdag_dso2  = gas.get(0, 1, 1);

                double d2lambda_dsw2  = d2lambdaw_dsw2  + d2lambdag_dsw2  + d2lambdao_dsw2;
                double d2lambda_dswso = d2lambdaw_dswso + d2lambdag_dswso + d2lambdao_dswso;
                double d2lambda_dsosw = d2lambda_dswso;
                double d2lambda_dso2  = d2lambdaw_dso2  + d2lambdag_dso2  + d2lambdao_dso2;

                double d2fw_dsw2  = inv_lambda*d2lambdaw_dsw2 - 2.0*inv_lambda2*dlambdaw_dsw*dlambda_dsw +
                                    twice_inv_lambda3*lambdaw*(dlambda_dsw*dlambda_dsw) - 
                                    lambdaw*inv_lambda2*d2lambda_dsw2;

                double d2fw_dswso = inv_lambda*d2lambdaw_dswso - inv_lambda2*dlambdaw_dsw*dlambda_dso -
                                    inv_lambda2*dlambdaw_dso*dlambda_dsw + 
                                    twice_inv_lambda3*lambdaw*dlambda_dso*dlambda_dsw -
                                    lambdaw*inv_lambda2*d2lambda_dswso;

                double d2fw_dsosw = d2fw_dswso;

                double d2fw_dso2  = inv_lambda*d2lambdaw_dso2 - 2.0*inv_lambda2*dlambdaw_dso*dlambda_dso +
                                    twice_inv_lambda3*lambdaw*(dlambda_dso*dlambda_dso) - 
                                    lambdaw*inv_lambda2*d2lambda_dso2;

                f.set(0, 0, 0, d2fw_dsw2);
                f.set(0, 0, 1, d2fw_dswso);
                f.set(0, 1, 0, d2fw_dsosw);
                f.set(0, 1, 1, d2fw_dso2);

                double d2fo_dsw2  = inv_lambda*d2lambdao_dsw2 - 2.0*inv_lambda2*dlambdao_dsw*dlambda_dsw +
                                    twice_inv_lambda3*lambdao*(dlambda_dsw*dlambda_dsw) - 
                                    lambdao*inv_lambda2*d2lambda_dsw2;

                double d2fo_dswso = inv_lambda*d2lambdao_dswso - inv_lambda2*dlambdao_dsw*dlambda_dso -
                                    inv_lambda2*dlambdao_dso*dlambda_dsw + 
                                    twice_inv_lambda3*lambdao*dlambda_dso*dlambda_dsw -
                                    lambdao*inv_lambda2*d2lambda_dswso;

                double d2fo_dsosw = d2fo_dswso;

                double d2fo_dso2  = inv_lambda*d2lambdao_dso2 - 2.0*inv_lambda2*dlambdao_dso*dlambda_dso +
                                    twice_inv_lambda3*lambdao*(dlambda_dso*dlambda_dso) - 
                                    lambdao*inv_lambda2*d2lambda_dso2;

                f.set(1, 0, 0, d2fo_dsw2);
                f.set(1, 0, 1, d2fo_dswso);
                f.set(1, 1, 0, d2fo_dsosw);
                f.set(1, 1, 1, d2fo_dso2);

            }
        }
    }

    // std::cout << "    Jet ends." << std::endl;

    return 2;
}


#include "ThreePhaseFlowImplicitHugoniotCurve.h"

ThreePhaseFlowImplicitHugoniotCurve::ThreePhaseFlowImplicitHugoniotCurve(ThreePhaseFlowSubPhysics *s): ImplicitHugoniotCurve(s->flux(), s->accumulation(), s->boundary()), 
                                                                                               subphysics_(s),
                                                                                               permeability_(s->permeability())  {
     method_ = IMPLICIT_HUGONIOT;
     info_ = std::string("ThreePhaseFlowImplicitHugoniotCurve");
}

ThreePhaseFlowImplicitHugoniotCurve::~ThreePhaseFlowImplicitHugoniotCurve(){
}

int ThreePhaseFlowImplicitHugoniotCurve::classical_function_on_square(ImplicitHugoniotCurve *obj, double *foncub, int i, int j){
    int is_square = obj->gridvalues()->cell_type(i, j);
    double dg1, dg2;
    double f_aux[4];
    double epsilon = 1.0e-5;

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            dg1 = obj->gridvalues()->G_on_grid(i + l, j + k).component(0) - obj->G_ref().component(0);
            dg2 = obj->gridvalues()->G_on_grid(i + l, j + k).component(1) - obj->G_ref().component(1);

            if (fabs(dg1) + fabs(dg2) >= epsilon) {
                double df1 = obj->gridvalues()->F_on_grid(i + l, j + k).component(0) - obj->F_ref().component(0);
                double df2 = obj->gridvalues()->F_on_grid(i + l, j + k).component(1) - obj->F_ref().component(1);

                f_aux[l * 2 + k] = dg2 * df1 - dg1*df2;
            } 
            else {
                // First-order expansion of F in terms of G.
                //

                double inv_det = 1.0 / (obj->JG_ref()(0) * obj->JG_ref()(3) - obj->JG_ref()(1) * obj->JG_ref()(2) );

                f_aux[l * 2 + k] = ((obj->JF_ref()(0) * obj->JG_ref()(3) - obj->JF_ref()(2) * obj->JG_ref()(1) + obj->JF_ref()(1) * obj->JG_ref()(2) - obj->JF_ref()(3) * obj->JG_ref()(0)) * dg1 * dg2 +
                        (obj->JF_ref()(1) * obj->JG_ref()(3) - obj->JF_ref()(3) * obj->JG_ref()(1)) * dg2 * dg2 +
                        (obj->JF_ref()(0) * obj->JG_ref()(2) - obj->JF_ref()(2) * obj->JG_ref()(0)) * dg1 * dg1) * inv_det;
            }
        }
    }

    foncub[1] = f_aux[0]; // Was: foncub[0][1]
    foncub[0] = f_aux[2]; // Was: foncub[0][0]
    foncub[3] = f_aux[1]; // Was: foncub[0][2]

    // Only useful if the cell is a square.
    //
    if (is_square == CELL_IS_SQUARE) foncub[2] = f_aux[3]; // Was: foncub[0][2]

    return 1;
}

int ThreePhaseFlowImplicitHugoniotCurve::threephase_function_on_square(ImplicitHugoniotCurve *obj, double *foncub, int i, int j) {
    int is_square = obj->gridvalues()->cell_type(i, j);

    RealVector state(2); //state(3);
    RealVector f_aux(4);

    double (*implicit_Hugoniot_function)(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state) = ((ThreePhaseFlowImplicitHugoniotCurve*)obj)->implicit_Hugoniot_function;

    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < 2; k++) {
            state(0) = obj->gridvalues()->grid(i + l, j + k).component(0);
            state(1) = obj->gridvalues()->grid(i + l, j + k).component(1);
/*            state(2) = 1.0 - state(0) - state(1);*/

            f_aux(l * 2 + k) = (*implicit_Hugoniot_function)((ThreePhaseFlowImplicitHugoniotCurve*)obj, state);
        }
    }

    foncub[1] = f_aux(0); // Was: foncub[0][1]
    foncub[0] = f_aux(2); // Was: foncub[0][0]
    foncub[3] = f_aux(1); // Was: foncub[0][2]

    // Only useful if the cell is a square.
    //
    if (is_square == CELL_IS_SQUARE) foncub[2] = f_aux(3); // Was: foncub[0][2]

    return 1;
}

int ThreePhaseFlowImplicitHugoniotCurve::function_on_square(double *foncub, int i, int j) {
//    std::cout << "Here 1" << std::endl;

    int info = (*fonsq)(this, foncub, i, j);

//    std::cout << "Here 2" << std::endl;

    return info;
}

double ThreePhaseFlowImplicitHugoniotCurve::gas_vertex_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state){
    RealVector redperm(3);
    obj->permeability_->reduced_permeability(state, redperm);

    JetMatrix water(2);
    obj->permeability_->PermeabilityWater_jet(state, 0, water);
    double lambda_w = water.get(0)/obj->muw; // kw/muw

    JetMatrix oil(2);
    obj->permeability_->PermeabilityOil_jet(state, 0, oil);
    double lambda_o = oil.get(0)/obj->muo; // ko/muo

    JetMatrix gas(2);
    obj->permeability_->PermeabilityGas_jet(state, 0, gas);

    JetMatrix mug_jet;
    obj->subphysics_->viscosity()->gas_viscosity_jet(state, 0, mug_jet);
    double mug = mug_jet.get(0);

    double lambda_g = gas.get(0)/mug; // kg/mug

    double oil_factor   = (redperm(1)/obj->muo)*(obj->vel + lambda_g*obj->rho_o_g + lambda_w*obj->rho_o_w);
    double water_factor = (redperm(0)/obj->muw)*(obj->vel + lambda_g*obj->rho_w_g + lambda_o*obj->rho_w_o);

    return oil_factor - water_factor;
}

double ThreePhaseFlowImplicitHugoniotCurve::water_vertex_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state){
    RealVector redperm(3);
    obj->permeability_->reduced_permeability(state, redperm);

    JetMatrix water(2);
    obj->permeability_->PermeabilityWater_jet(state, 0, water);
    double lambda_w = water.get(0)/obj->muw; // kw/muw

    JetMatrix oil(2);
    obj->permeability_->PermeabilityOil_jet(state, 0, oil);
    double lambda_o = oil.get(0)/obj->muo; // ko/muo

    JetMatrix gas(2);
    obj->permeability_->PermeabilityGas_jet(state, 0, gas);

    JetMatrix mug_jet;
    obj->subphysics_->viscosity()->gas_viscosity_jet(state, 0, mug_jet);
    double mug = mug_jet.get(0);

    double lambda_g = gas.get(0)/mug; // kg/mug

    double oil_factor   = (redperm(1)/obj->muo)*(obj->vel + lambda_g*obj->rho_o_g + lambda_w*obj->rho_o_w);
    double gas_factor   = (redperm(2)/mug)*(obj->vel + lambda_w*obj->rho_g_w + lambda_o*obj->rho_g_o);

    return oil_factor - gas_factor;
}

double ThreePhaseFlowImplicitHugoniotCurve::oil_vertex_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state){
    RealVector redperm(3);
    obj->permeability_->reduced_permeability(state, redperm);

    JetMatrix water(2);
    obj->permeability_->PermeabilityWater_jet(state, 0, water);
    double lambda_w = water.get(0)/obj->muw; // kw/muw

    JetMatrix oil(2);
    obj->permeability_->PermeabilityOil_jet(state, 0, oil);
    double lambda_o = oil.get(0)/obj->muo; // ko/muo

    JetMatrix gas(2);
    obj->permeability_->PermeabilityGas_jet(state, 0, gas);

    JetMatrix mug_jet;
    obj->subphysics_->viscosity()->gas_viscosity_jet(state, 0, mug_jet);
    double mug = mug_jet.get(0);

    double lambda_g = gas.get(0)/mug; // kg/mug

    double water_factor = (redperm(0)/obj->muw)*(obj->vel + lambda_g*obj->rho_w_g + lambda_o*obj->rho_w_o);
    double gas_factor   = (redperm(2)/mug)*(obj->vel + lambda_w*obj->rho_g_w + lambda_o*obj->rho_g_o);

    return water_factor - gas_factor;
}

double ThreePhaseFlowImplicitHugoniotCurve::water_gas_side_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state){
    double sw = state(0);
    double sw_ref = obj->reference_point.point(0);

    RealVector redperm(3);
    obj->permeability_->reduced_permeability(state, redperm);

    JetMatrix water(2);
    obj->permeability_->PermeabilityWater_jet(state, 0, water);
    double lambda_w = water.get(0)/obj->muw; // kw/muw

    JetMatrix oil(2);
    obj->permeability_->PermeabilityOil_jet(state, 0, oil);
    double lambda_o = oil.get(0)/obj->muo; // ko/muo

    JetMatrix gas(2);
    obj->permeability_->PermeabilityGas_jet(state, 0, gas);

    JetMatrix mug_jet;
    obj->subphysics_->viscosity()->gas_viscosity_jet(state, 0, mug_jet);
    double mug = mug_jet.get(0);

    double lambda_g = gas.get(0)/mug; // kg/mug

    double lambda = lambda_w + lambda_o + lambda_g;

    double term1   = (redperm(1)/(obj->muo*lambda))*(obj->vel + lambda_g*obj->rho_o_g + lambda_w*obj->rho_o_w)*(sw_ref - sw);

    double term2   = (obj->lambda_w_ref/(obj->lambda_w_ref + obj->lambda_g_ref))*(obj->vel + obj->lambda_g_ref*obj->rho_w_g);

    double term3 = (lambda_w/lambda)*(obj->vel + lambda_g*obj->rho_w_g + lambda_o*obj->rho_w_o);

    return term1 - term2 + term3;
}

double ThreePhaseFlowImplicitHugoniotCurve::water_oil_side_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state){
    double sw = state(0);
    double sw_ref = obj->reference_point.point(0);

    RealVector redperm(3);
    obj->permeability_->reduced_permeability(state, redperm);

    JetMatrix water(2);
    obj->permeability_->PermeabilityWater_jet(state, 0, water);
    double lambda_w = water.get(0)/obj->muw; // kw/muw

    JetMatrix oil(2);
    obj->permeability_->PermeabilityOil_jet(state, 0, oil);
    double lambda_o = oil.get(0)/obj->muo; // ko/muo

    JetMatrix gas(2);
    obj->permeability_->PermeabilityGas_jet(state, 0, gas);

    JetMatrix mug_jet;
    obj->subphysics_->viscosity()->gas_viscosity_jet(state, 0, mug_jet);
    double mug = mug_jet.get(0);

    double lambda_g = gas.get(0)/mug; // kg/mug

    double lambda = lambda_w + lambda_o + lambda_g;

    double term1   = (redperm(2)/(mug*lambda))*(obj->vel + lambda_w*obj->rho_g_w + lambda_o*obj->rho_g_o)*(sw_ref - sw);

    double term2 = (obj->lambda_w_ref/(obj->lambda_w_ref + obj->lambda_o_ref))*(obj->vel + obj->lambda_o_ref*obj->rho_w_o);

    double term3   = (lambda_w/lambda)*(obj->vel + lambda_g*obj->rho_w_g + lambda_o*obj->rho_w_o);

    return term1 - term2 + term3;
}

double ThreePhaseFlowImplicitHugoniotCurve::oil_gas_side_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state){
    double so = state(1);
    double so_ref = obj->reference_point.point(1);

    RealVector redperm(3);
    obj->permeability_->reduced_permeability(state, redperm);

    JetMatrix water(2);
    obj->permeability_->PermeabilityWater_jet(state, 0, water);
    double lambda_w = water.get(0)/obj->muw; // kw/muw

    JetMatrix oil(2);
    obj->permeability_->PermeabilityOil_jet(state, 0, oil);
    double lambda_o = oil.get(0)/obj->muo; // ko/muo

    JetMatrix gas(2);
    obj->permeability_->PermeabilityGas_jet(state, 0, gas);

    JetMatrix mug_jet;
    obj->subphysics_->viscosity()->gas_viscosity_jet(state, 0, mug_jet);
    double mug = mug_jet.get(0);

    double lambda_g = gas.get(0)/mug; // kg/mug

    double lambda = lambda_w + lambda_o + lambda_g;

    double term1 = (redperm(0)/(obj->muw*lambda))*(obj->vel + lambda_g*obj->rho_w_g + lambda_o*obj->rho_w_o)*(so_ref - so);

    double term2 = (obj->lambda_o_ref/(obj->lambda_o_ref + obj->lambda_g_ref))*(obj->vel + obj->lambda_g_ref*obj->rho_o_g);

    double term3 = (lambda_o/lambda)*(obj->vel + lambda_g*obj->rho_o_g + lambda_w*obj->rho_o_w);

    return term1 - term2 + term3;
}

void ThreePhaseFlowImplicitHugoniotCurve::curve(const ReferencePoint &ref, int type, std::vector<Curve> &c){
    if (subphysics_ != 0) gv = subphysics_->gridvalues();

    if (type == THREEPHASEFLOWIMPLICITHUGONIOTCURVE_G_VERTEX){
        implicit_Hugoniot_function = &gas_vertex_function;
        fonsq = &threephase_function_on_square;

        // No projection is necessary.
        reference_point.point = subphysics_->G();

        // Add GW
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->G(), subphysics_->W(), n, temp_curve);
            c.push_back(temp_curve);
        }

        // Add GO
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->G(), subphysics_->O(), n, temp_curve);
            c.push_back(temp_curve);
        }

//                // Add GW side as two Curves.
//                //
//                {
//                    Curve temp_curve;
//                    Utilities::regularly_sampled_segment(subphysics_->G(), point_on_side[1], n, temp_curve);
//                    c.push_back(temp_curve);

//                    Utilities::regularly_sampled_segment(point_on_side[1], subphysics_->W(), n, temp_curve);
//                    c.push_back(temp_curve);
//                }

//                // Add GO side as two Curves.
//                //
//                {
//                    Curve temp_curve;
//                    Utilities::regularly_sampled_segment(subphysics_->G(), point_on_side[0], n, temp_curve);
//                    c.push_back(temp_curve);

//                    Utilities::regularly_sampled_segment(point_on_side[0], subphysics_->O(), n, temp_curve);
//                    c.push_back(temp_curve);
//                }
    }
    else if (type == THREEPHASEFLOWIMPLICITHUGONIOTCURVE_W_VERTEX){
        implicit_Hugoniot_function = &water_vertex_function;
        fonsq = &threephase_function_on_square;

        // No projection is necessary.
        reference_point.point = subphysics_->W();

        // Add GW
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->G(), subphysics_->W(), n, temp_curve);
            c.push_back(temp_curve);
        }

        // Add WO
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->W(), subphysics_->O(), n, temp_curve);
            c.push_back(temp_curve);
        }

//                // Add GW side as two Curves.
//                //
//                {
//                    Curve temp_curve;
//                    Utilities::regularly_sampled_segment(subphysics_->G(), point_on_side[1], n, temp_curve);
//                    c.push_back(temp_curve);

//                    Utilities::regularly_sampled_segment(point_on_side[1], subphysics_->W(), n, temp_curve);
//                    c.push_back(temp_curve);
//                }

//                // Add WO side as two Curves.
//                //
//                {
//                    Curve temp_curve;
//                    Utilities::regularly_sampled_segment(subphysics_->W(), point_on_side[2], n, temp_curve);
//                    c.push_back(temp_curve);

//                    Utilities::regularly_sampled_segment(point_on_side[2], subphysics_->O(), n, temp_curve);
//                    c.push_back(temp_curve);
//                }
    }
    else if (type == THREEPHASEFLOWIMPLICITHUGONIOTCURVE_O_VERTEX){
        implicit_Hugoniot_function = &oil_vertex_function;
        fonsq = &threephase_function_on_square;

        // No projection is necessary.
        reference_point.point = subphysics_->O();

        // Add GO
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->G(), subphysics_->O(), n, temp_curve);
            c.push_back(temp_curve);
        }

        // Add WO
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->W(), subphysics_->O(), n, temp_curve);
            c.push_back(temp_curve);
        }

//                // Add GO side as two Curves.
//                //
//                {
//                    Curve temp_curve;
//                    Utilities::regularly_sampled_segment(subphysics_->G(), point_on_side[0], n, temp_curve);
//                    c.push_back(temp_curve);

//                    Utilities::regularly_sampled_segment(point_on_side[0], subphysics_->O(), n, temp_curve);
//                    c.push_back(temp_curve);
//                }

//                // Add WO side as two Curves.
//                //
//                {
//                    Curve temp_curve;
//                    Utilities::regularly_sampled_segment(subphysics_->W(), point_on_side[2], n, temp_curve);
//                    c.push_back(temp_curve);

//                    Utilities::regularly_sampled_segment(point_on_side[2], subphysics_->O(), n, temp_curve);
//                    c.push_back(temp_curve);
//                }
    }
    else if (type == THREEPHASEFLOWIMPLICITHUGONIOTCURVE_GW_SIDE) {
        implicit_Hugoniot_function = &water_gas_side_function;
        fonsq = &threephase_function_on_square;

        RealVector projected_point = project_point_onto_line_2D(ref.point, subphysics_->G(), subphysics_->W());
        reference_point = ReferencePoint(projected_point, subphysics_->flux(), subphysics_->accumulation(), 0);

        // Add GW
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->G(), subphysics_->W(), n, temp_curve);
            c.push_back(temp_curve);
        }
    }
    else if (type == THREEPHASEFLOWIMPLICITHUGONIOTCURVE_WO_SIDE) {
        implicit_Hugoniot_function = &water_oil_side_function;
        fonsq = &threephase_function_on_square;

        RealVector projected_point = project_point_onto_line_2D(ref.point, subphysics_->W(), subphysics_->O());
        reference_point = ReferencePoint(projected_point, subphysics_->flux(), subphysics_->accumulation(), 0);

        // Add WO
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->W(), subphysics_->O(), n, temp_curve);
            c.push_back(temp_curve);
        }
    }
    else if (type == THREEPHASEFLOWIMPLICITHUGONIOTCURVE_GO_SIDE) {
        implicit_Hugoniot_function = &oil_gas_side_function;
        fonsq = &threephase_function_on_square;

        RealVector projected_point = project_point_onto_line_2D(ref.point, subphysics_->G(), subphysics_->O());
        reference_point = ReferencePoint(projected_point, subphysics_->flux(), subphysics_->accumulation(), 0);

        // Add GO
        //
        {
            int n = 100;
            Curve temp_curve;
            Utilities::regularly_sampled_segment(subphysics_->G(), subphysics_->O(), n, temp_curve);
            c.push_back(temp_curve);
        }
    }
    else {
        // Generic.
        //
        fonsq = &classical_function_on_square;

        reference_point = ref;

        // Invoke the base's class' method.
        //
        ImplicitHugoniotCurve::curve(ref, IMPLICITHUGONIOTCURVE_GENERIC_POINT, c);

        return;
    }

    vel = subphysics_->vel()->value();
    muw = subphysics_->muw()->value();
    muo = subphysics_->muo()->value();
//    mug = subphysics_->mug()->value();
    grw = subphysics_->grw()->value();
    gro = subphysics_->gro()->value();
    grg = subphysics_->grg()->value();

    rho_w_o = grw - gro; // rho12
    rho_o_w = -rho_w_o;  // rho21

    rho_w_g = grw - grg; // rho13
    rho_g_w = -rho_w_g;  // rho31

    rho_o_g = gro - grg; // rho23
    rho_g_o = -rho_o_g;  // rho32

    JetMatrix water(2);
    permeability_->PermeabilityWater_jet(reference_point.point, 0, water);
    lambda_w_ref = water.get(0)/muw; // kw/muw

    JetMatrix oil(2);
    permeability_->PermeabilityOil_jet(reference_point.point, 0, oil);
    lambda_o_ref = oil.get(0)/muo; // ko/muo

    JetMatrix gas(2);
    permeability_->PermeabilityGas_jet(reference_point.point, 0, gas);

    JetMatrix mug_jet;
    subphysics_->viscosity()->gas_viscosity_jet(reference_point.point, 0, mug_jet);
    double mug = mug_jet.get(0);

    lambda_g_ref = gas.get(0)/mug; // kg/mug

    std::vector<RealVector> hugoniot_curve;
    std::vector< std::deque <RealVector> > curves;
    std::vector <bool> is_circular;

    int method = SEGMENTATION_METHOD;
    int info = ContourMethod::contour2d(this, hugoniot_curve, curves, is_circular, method);

    for (int i = 0; i < hugoniot_curve.size()/2; i++){
        Curve temp;
        temp.curve.push_back(hugoniot_curve[2*i]);
        temp.curve.push_back(hugoniot_curve[2*i + 1]);

        c.push_back(temp);
    }

    return;
}


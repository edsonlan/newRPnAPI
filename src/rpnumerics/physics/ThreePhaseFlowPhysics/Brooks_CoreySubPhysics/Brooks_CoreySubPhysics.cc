#include "Brooks_CoreySubPhysics.h"

Brooks_CoreySubPhysics::Brooks_CoreySubPhysics() : ThreePhaseFlowSubPhysics(){
    // Every ThreePhaseFlowSubPhysics has these.
    //
    muw_parameter = new Parameter(std::string("muw"), 1.0);
    muo_parameter = new Parameter(std::string("muo"), 2.0);
    mug_parameter = new Parameter(std::string("mug"), 0.5);

    // Brooks_CoreySubPhysics-specific.
    //
    lambda_parameter = new Parameter(std::string("lambda"), .2, 7.0, 1.0);
    cnw_parameter = new Parameter(std::string("cnw"), 0.0);
    cno_parameter = new Parameter(std::string("cno"), 0.0);
    cng_parameter = new Parameter(std::string("cng"), 0.0);

    // Not in the original Brooks-Corey, but needed by the permeability.
    // TODO: Decide if they are to be added to the list of
    //       parameters. They could be NOT added, so they
    //       would not be modified.
    //
    grw_parameter = new Parameter(std::string("grw"), 1.0);
    gro_parameter = new Parameter(std::string("gro"), 1.0);
    grg_parameter = new Parameter(std::string("grg"), 1.0);
    vel_parameter = new Parameter(std::string("vel"), 1.0);

    equation_parameter_.push_back(muw_parameter);
    equation_parameter_.push_back(muo_parameter);
    equation_parameter_.push_back(mug_parameter);
    
    
    
    
    
    

    equation_parameter_.push_back(grw_parameter); // TODO: Add or not?
    equation_parameter_.push_back(gro_parameter); // TODO: Add or not?
    equation_parameter_.push_back(grg_parameter); // TODO: Add or not?
    equation_parameter_.push_back(vel_parameter); // TODO: Add or not?

    // Permeability.
    //
    permeability_ = new Brooks_CoreyPermeability(lambda_parameter, cnw_parameter, cno_parameter, cng_parameter, this);
    auxiliaryfunctions_.push_back(permeability_);

    // Flux.
    //
    flux_ = new Brooks_CoreyFluxFunction(muw_parameter, muo_parameter, mug_parameter, (Brooks_CoreyPermeability*)permeability_);

    // GridValues.
    //
    std::vector<int> number_of_cells(2);
    number_of_cells[0] = 128;
    number_of_cells[1] = 128;

    gridvalues_ = new GridValues(boundary_, boundary_->minimums(), boundary_->maximums(), number_of_cells);

    // HugoniotContinuation.
    //
//    hugoniotcontinuation_ = new HugoniotContinuation2D2D(flux_, accumulation_, boundary_);

    // ImplicitHugoniot.
    //
    ImplicitHugoniotCurve *ihc = new ImplicitHugoniotCurve(flux_, accumulation_, boundary_);
    ihc->subphysics(this);

    hugoniot_curve.push_back(ihc);

    // ThreePhaseFlowImplicitHugoniot.
    //
    ThreePhaseFlowImplicitHugoniotCurve *tpfih = new ThreePhaseFlowImplicitHugoniotCurve(this);
    hugoniot_curve.push_back(tpfih);

    // Rarefaction.
    //
    rarefactioncurve_ = new RarefactionCurve(accumulation_, flux_, boundary_);

    // Inflection.
    //
    inflection_curve_ = new Inflection_Curve;

    // HugoniotContinuation.
    //
    hugoniotcontinuation_ = new HugoniotContinuation2D2D(flux_, accumulation_, boundary_);
//    hugoniotcontinuation_ = new HugoniotContinuation_nDnD(flux_, accumulation_, boundary_);
    
    hugoniot_curve.push_back(hugoniotcontinuation_);

    shockcurve_ = new ShockCurve(hugoniotcontinuation_);

    // ODE solver.
    //
    odesolver_ = new LSODE;

    // Composite.
    //
    compositecurve_ = new CompositeCurve(accumulation_, flux_, boundary_, shockcurve_, 0);

    // WaveCurve.
    //
    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);

    // Viscosity.
    //
    viscosity_ = new Brooks_CoreyViscosity(this);

    // Mobility.
    //
    mobility_ = new ThreePhaseFlowMobility(this);

    // Info.
    //
    info_subphysics_ = std::string("Brooks-Corey");
}

Brooks_CoreySubPhysics::~Brooks_CoreySubPhysics(){
    delete mobility_;

    delete viscosity_;

    delete wavecurvefactory_;

    delete compositecurve_;

    delete odesolver_;

    delete shockcurve_;

    delete inflection_curve_;

    delete rarefactioncurve_;

    for (int i = 0; i < hugoniot_curve.size(); i++) delete hugoniot_curve[i];

    delete gridvalues_;

    delete flux_;

    for (int i = 0; i < equation_parameter_.size(); i++) delete equation_parameter_[i];

    for (int i = 0; i < auxiliaryfunctions_.size(); i++) delete auxiliaryfunctions_[i];
}

void Brooks_CoreySubPhysics::shock_cases(std::vector<int> &type, std::vector<std::string> &name) const {
    type.clear();
    name.clear();

    type.push_back(BROOKS_COREYGENERICPOINT);
    name.push_back(std::string("Generic point"));

    return;
}


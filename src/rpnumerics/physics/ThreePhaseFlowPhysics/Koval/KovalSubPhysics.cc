#include "KovalSubPhysics.h"

KovalSubPhysics::KovalSubPhysics() : ThreePhaseFlowSubPhysics(){
    grw_parameter = new Parameter(std::string("grw"), 1.0);
    gro_parameter = new Parameter(std::string("gro"), 1.0);
    grg_parameter = new Parameter(std::string("grg"), 1.0);

    muw_parameter = new Parameter(std::string("muw"), 1.0 + 1e-3);
    muo_parameter = new Parameter(std::string("muo"), 1.0 + 2e-3);
    mug_parameter = new Parameter(std::string("mug"), 1.0 - 3e-3);

    vel_parameter = new Parameter(std::string("vel"), 1.0);

    equation_parameter_.push_back(grw_parameter);
    equation_parameter_.push_back(gro_parameter);
    equation_parameter_.push_back(grg_parameter);

    equation_parameter_.push_back(muw_parameter);
    equation_parameter_.push_back(muo_parameter);
    equation_parameter_.push_back(mug_parameter);

    equation_parameter_.push_back(vel_parameter);

    // Permeability.
    //
    permeability_ = new KovalPermeability(this);

    // Flux.
    //
    flux_ = new KovalFluxFunction(grw_parameter, gro_parameter, grg_parameter,
                                  muw_parameter, muo_parameter, mug_parameter,
                                  vel_parameter);

    // GridValues.
    //
    std::vector<int> number_of_cells(2);
    number_of_cells[0] = 128;
    number_of_cells[1] = 128;

    gridvalues_ = new GridValues(boundary_, boundary_->minimums(), boundary_->maximums(), number_of_cells);

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

    // Shock curve.
    //
    hugoniotcontinuation_ = new HugoniotContinuation2D2D(flux_, accumulation_, boundary_);
//    hugoniot_curve.push_back(hugoniotcontinuation_);
    shockcurve_ = new ShockCurve(hugoniotcontinuation_);

    // Composite.
//    Stone_Explicit_Bifurcation_Curves bc((StoneFluxFunction*)flux);
    compositecurve_ = new CompositeCurve(accumulation_, flux_, boundary_, shockcurve_, 0/*&bc*/);

    odesolver_ = new LSODE;

    // WaveCurve.
    //
    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);

    // Inflection.
    //
    inflection_curve_ = new Inflection_Curve;

    // Viscosity.
    //
    viscosity_ = new KovalViscosity(this);

    // Mobility.
    //
    mobility_ = new ThreePhaseFlowMobility(this);

    info_subphysics_ = std::string("Koval");
}

KovalSubPhysics::~KovalSubPhysics(){
    delete mobility_;
    delete viscosity_;
    delete inflection_curve_;
    delete wavecurvefactory_;
    delete odesolver_;
    delete compositecurve_;
    delete shockcurve_;
    delete hugoniotcontinuation_;
    delete rarefactioncurve_;
    for (int i = 0; i < hugoniot_curve.size(); i++) delete hugoniot_curve[i];

    for (int i = 0; i < equation_parameter_.size(); i++) delete equation_parameter_[i];

    delete flux_;

    delete permeability_;
}


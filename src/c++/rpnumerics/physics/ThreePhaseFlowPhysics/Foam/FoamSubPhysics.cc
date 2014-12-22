#include "FoamSubPhysics.h"

FoamSubPhysics::FoamSubPhysics() : ThreePhaseFlowSubPhysics(){
    muw_parameter = new Parameter(std::string("muw"), 1.0);
    muo_parameter = new Parameter(std::string("muo"), 1.0);
    mug_parameter = new Parameter(std::string("mug0"), 1.0); // Notice that, even though the name is changed, this parameter is just mug.

    grw_parameter = new Parameter(std::string("grw"), 1.0);
    gro_parameter = new Parameter(std::string("gro"), 1.0);
    grg_parameter = new Parameter(std::string("grg"), 1.0);

    vel_parameter = new Parameter(std::string("vel"), 1.0);

    cnw_parameter = new Parameter(std::string("cnw"), 0.0);
    cno_parameter = new Parameter(std::string("cno"), 0.0);
    cng_parameter = new Parameter(std::string("cng"), 0.0);

    nw_parameter = new Parameter(std::string("nw"), 2.0);
    no_parameter = new Parameter(std::string("no"), 2.0);
    ng_parameter = new Parameter(std::string("ng"), 2.0);

    equation_parameter_.push_back(muw_parameter);
    equation_parameter_.push_back(muo_parameter);
    equation_parameter_.push_back(mug_parameter);

    equation_parameter_.push_back(grw_parameter);
    equation_parameter_.push_back(gro_parameter);
    equation_parameter_.push_back(grg_parameter);

    equation_parameter_.push_back(vel_parameter);

    equation_parameter_.push_back(cnw_parameter);
    equation_parameter_.push_back(cno_parameter);
    equation_parameter_.push_back(cng_parameter);

    equation_parameter_.push_back(nw_parameter);
    equation_parameter_.push_back(no_parameter);
    equation_parameter_.push_back(ng_parameter);

    // TODO: These values must be set correctly!
    epdry = new Parameter(std::string("epdry"), 0.0);
    fdry  = new Parameter(std::string("fdry"), 0.0);
    foil  = new Parameter(std::string("foil"), 0.0);
    fmdry = new Parameter(std::string("fmdry"), 0.0);
    fmmob = new Parameter(std::string("fmmob"), 55000.0);
    fmoil = new Parameter(std::string("fmoil"), 0.0);

    equation_parameter_.push_back(epdry);
    equation_parameter_.push_back(fdry);
    equation_parameter_.push_back(foil);
    equation_parameter_.push_back(fmdry);
    equation_parameter_.push_back(fmmob);
    equation_parameter_.push_back(fmoil);

    // Permeability.
    //
    permeability_ = new FoamPermeability(cnw_parameter, cno_parameter, cng_parameter,
                                         nw_parameter, no_parameter, ng_parameter,
                                         this);

    // Gas viscosity.
    //
    viscosity_ = new FoamViscosity(mug_parameter, 
                                   epdry,
                                   fdry,
                                   foil,
                                   fmdry,
                                   fmmob,
                                   fmoil,
                                   this);

    // Flux.
    //
    flux_ = new FoamFluxFunction(grw_parameter, gro_parameter, grg_parameter,
                                 muw_parameter, muo_parameter,
                                 vel_parameter,
                                 (FoamPermeability*)permeability_,
                                 (FoamViscosity*)viscosity_);

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
    hugoniot_curve.push_back(hugoniotcontinuation_);

    shockcurve_ = new ShockCurve(hugoniotcontinuation_);

    // Composite.
//    Stone_Explicit_Bifurcation_Curves bc((StoneFluxFunction*)flux);
    compositecurve_ = new CompositeCurve(accumulation_, flux_, boundary_, shockcurve_, 0/*&bc*/);

    odesolver_ = new LSODE;

    // WaveCurve.
    //
//    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);
    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);

    // Inflection.
    //
    inflection_curve_ = new Inflection_Curve;

//    // BifurcationCurve.
//    //
//    bifurcationcurve_ = new FoamTransitionalLine(this);

    // Mobility.
    //
    mobility_ = new ThreePhaseFlowMobility(this);

    info_subphysics_ = std::string("Foam");
}

FoamSubPhysics::~FoamSubPhysics(){
    delete mobility_;
    delete inflection_curve_;
    delete wavecurvefactory_;
    delete odesolver_;
    delete compositecurve_;
    delete shockcurve_;
    delete hugoniotcontinuation_;
    delete rarefactioncurve_;
    for (int i = 0; i < hugoniot_curve.size(); i++) delete hugoniot_curve[i];

    delete flux_;

    delete viscosity_;

    delete permeability_;

    for (int i = 0; i < equation_parameter_.size(); i++) delete equation_parameter_[i];
}


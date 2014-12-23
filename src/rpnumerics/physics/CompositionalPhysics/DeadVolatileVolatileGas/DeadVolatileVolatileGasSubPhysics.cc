#include "DeadVolatileVolatileGasSubPhysics.h"

DeadVolatileVolatileGasSubPhysics::DeadVolatileVolatileGasSubPhysics() : SubPhysics() {
    // Equation parameters.
    //
    re_parameter  = new Parameter(std::string("re"), 44.0);// re_parameter->value(0.0);
    rg_parameter  = new Parameter(std::string("rg"), 100.0);
    phi_parameter = new Parameter(std::string("phi"), 0.2);
    
    // Thermodynamics parameters.
    //
    B_parameter = new Parameter(std::string("B"), 7.0);
    D_parameter = new Parameter(std::string("D"), 3.0);

    mu_oB_parameter = new Parameter(std::string("mu_oB"), 10.0);
    mu_oD_parameter = new Parameter(std::string("mu_oD"), 8.0);
    mu_G_parameter  = new Parameter(std::string("mu_G"), 1.0);

    // Thermodynamics proper.
    //
    thermo = new DeadVolatileVolatileGasThermodynamics(B_parameter, D_parameter, mu_oB_parameter, mu_oD_parameter, mu_G_parameter, rg_parameter, re_parameter);
    auxiliaryfunctions_.push_back(thermo);
    
    // Hydrodynamics proper.
    //
    hydro  = new DeadVolatileVolatileGasHydrodynamics(thermo);
    auxiliaryfunctions_.push_back(hydro);

//    equation_parameter_.push_back(re_parameter);
//    equation_parameter_.push_back(rg_parameter);
    equation_parameter_.push_back(phi_parameter);

    // Flux, accumulation.
    //
    flux_         = new DeadVolatileVolatileGasFluxFunction(thermo, hydro);
    accumulation_ = new DeadVolatileVolatileGasAccumulationFunction(phi_parameter, thermo);

    // Boundary.
    //
    double eps = 1e-5;
    RealVector min(3), max(3);
    min(0) = 0.0 + eps;
    min(1) = 0.0 + eps; 
    min(2) = -1.0;

    max(0) = 1.0 - eps;
    max(1) = 1.0 - eps;
    max(2) = 1.0;

    boundary_     = new RectBoundary(min, max);

    // GridValues.
    //
    std::vector<int> number_of_cells(2);
    number_of_cells[0] = 128;
    number_of_cells[1] = 128;

    gridvalues_ = new GridValues(boundary_, boundary_->minimums(), boundary_->maximums(), number_of_cells);

    // Implicit Hugoniot.
    //
    Hugoniot_TP *ihc = new Hugoniot_TP(flux_, accumulation_, boundary_);
    ihc->subphysics(this);

    hugoniot_curve.push_back(ihc);

    // Hugoniot continuation

//    hugoniotcontinuation_ = new HugoniotContinuation3D2D(flux_, accumulation_, boundary_);
    hugoniotcontinuation_ = new HugoniotContinuation_nDnD(flux_, accumulation_, boundary_);

    hugoniot_curve.push_back(hugoniotcontinuation_);

    // Coincidence, Coincidence contour.
    //
    coincidence_ = new DeadVolatileVolatileGasCoincidence(thermo, hydro,
                                                     re_parameter, rg_parameter, 
                                                     phi_parameter);

    coincidence_contour_ = new Coincidence_Contour(coincidence_);

    // Rarefaction.
    //
    rarefactioncurve_ = new RarefactionCurve(accumulation_, flux_, boundary_);

    // Inflection.
    //
    inflection_curve_ = new Inflection_Curve;

    evap_ = new DeadVolatileVolatileGasEvaporationExtension((const DeadVolatileVolatileGasFluxFunction*)flux_, 
                                                    (const DeadVolatileVolatileGasAccumulationFunction*)accumulation_, 
                                                    (DeadVolatileVolatileGasCoincidence*)coincidence_,
                                                    phi_parameter);
    extension_curve.push_back(evap_);

    // Shockcurve.
    //
    shockcurve_ = new ShockCurve(hugoniotcontinuation_);

    // Composite.
    //
    compositecurve_ = new DeadVolatileVolatileGasCompositeCurve(evap_, accumulation_, flux_, boundary_, shockcurve_, 0/*&bc*/);

    odesolver_ = new LSODE;

    // Wavecurvefactory.
    //
    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);

    // Canvas-related.
    //
    transformation_matrix_ = DoubleMatrix::eye(2);

    xlabel_ = std::string("s");
    ylabel_ = std::string("y");
    
    // Info.
    //
    info_subphysics_ = std::string("DeadVolatileVolatileGasSubPhysics");
}

DeadVolatileVolatileGasSubPhysics::~DeadVolatileVolatileGasSubPhysics(){
    delete wavecurvefactory_;

    delete odesolver_;

    delete compositecurve_;

    delete shockcurve_;

    for (int i = 0; i < extension_curve.size(); i++) delete extension_curve[i];

    delete inflection_curve_;

    delete rarefactioncurve_;

    delete coincidence_contour_;
    delete coincidence_;
//    delete godcoincidence_;

    // Not needed, will be deleted by the block below.
//    delete hugoniotcontinuation_;

    // Not sure if this should really be done like this.
    // Perhaps it is best to eliminate only the HugoniotCurves that were instantiated
    // in this class, and let the rest be deleted at the father's dtor. 
    //
    for (int i = 0; i < hugoniot_curve.size(); i++) delete hugoniot_curve[i];

    delete gridvalues_;

    delete boundary_;
    delete accumulation_;
    delete flux_;

    delete rg_parameter;
    delete re_parameter;

    delete hydro;
    delete thermo;

    delete mu_G_parameter;
    delete mu_oD_parameter;
    delete mu_oB_parameter;

    delete D_parameter;
    delete B_parameter;
}

void DeadVolatileVolatileGasSubPhysics::shock_cases(std::vector<int> &type, std::vector<std::string> &name) const {
    type.clear();
    name.clear();

    type.push_back(DEADVOLATILEVOLATILEGAS_GENERIC_POINT);
    name.push_back(std::string("Generic point"));

    return;
}

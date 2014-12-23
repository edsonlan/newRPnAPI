#include "ICDOWSubPhysics.h"

ICDOWSubPhysics::ICDOWSubPhysics() : SubPhysics() {
    // Equation parameters.
    //
    muw_parameter = new Parameter(std::string("muw"), 1e-3); // Pa*s
    muo_parameter = new Parameter(std::string("muo"), 2e-3); // Pa*s

    grw_parameter = new Parameter(std::string("grw"), 0.0);
    gro_parameter = new Parameter(std::string("gro"), 0.0);

    vel_parameter = new Parameter(std::string("vel"), 1.0);

    swc_parameter = new Parameter(std::string("swc"), 0.15); // m^3/m^3

    lambda_parameter = new Parameter(std::string("lambda"), 2.0);

    phi_parameter = new Parameter(std::string("phi"), .2); // m^3/m^3

    equation_parameter_.push_back(muw_parameter);
    equation_parameter_.push_back(muo_parameter);
    equation_parameter_.push_back(grw_parameter);
    equation_parameter_.push_back(gro_parameter);
    equation_parameter_.push_back(vel_parameter);
    equation_parameter_.push_back(swc_parameter);
    equation_parameter_.push_back(lambda_parameter);
    equation_parameter_.push_back(phi_parameter);

    // Hydrodynamics.
    //
    hydrodynamics_  = new ICDOWHydrodynamics(muw_parameter, muo_parameter,
                                             grw_parameter, gro_parameter,
                                             vel_parameter,
                                             swc_parameter, lambda_parameter);

//    auxiliaryfunctions_.push_back(hydro);

    // Chemistry.
    //
    chemistry_ = new ICDOWChemistry();

//    auxiliaryfunctions_.push_back(chemistry);

    // Flux, accumulation.
    //
    flux_         = new ICDOWFluxFunction(chemistry_, hydrodynamics_);
    accumulation_ = new ICDOWAccumulationFunction(phi_parameter, chemistry_);

    // Boundary.
    //
    double eps = 1e-5;
    RealVector min(3), max(3);
    min(0) = 0.0 + eps;
    min(1) = -15.0; 
    min(2) = -1.0;

    max(0) = 1.0 - eps;
    max(1) = -7.0;
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

//    // Hugoniot continuation

////    hugoniotcontinuation_ = new HugoniotContinuation3D2D(flux_, accumulation_, boundary_);
//    hugoniotcontinuation_ = new HugoniotContinuation_nDnD(flux_, accumulation_, boundary_);

//    hugoniot_curve.push_back(hugoniotcontinuation_);

    // Coincidence, Coincidence contour.
    //
    coincidence_ = new ICDOWCoincidence(hydrodynamics_, phi_parameter);

    coincidence_contour_ = new Coincidence_Contour(coincidence_);

    // Rarefaction.
    //
    rarefactioncurve_ = new RarefactionCurve(accumulation_, flux_, boundary_);

    // Inflection.
    //
    inflection_curve_ = new Inflection_Curve;

//    evap_ = new ICDOWEvaporationExtension((const ICDOWFluxFunction*)flux_, 
//                                                    (const ICDOWAccumulationFunction*)accumulation_, 
//                                                    (ICDOWCoincidence*)coincidence_,
//                                                    phi_parameter);
//    extension_curve.push_back(evap_);

//    // Shockcurve.
//    //
//    shockcurve_ = new ShockCurve(hugoniotcontinuation_);

//    // Composite.
//    //
//    compositecurve_ = new ICDOWCompositeCurve(evap_, accumulation_, flux_, boundary_, shockcurve_, 0/*&bc*/);

//    odesolver_ = new LSODE;

//    // Wavecurvefactory.
//    //
//    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);

    // Canvas-related.
    //
    transformation_matrix_ = DoubleMatrix::eye(2);

    xlabel_ = std::string("sw");
    ylabel_ = std::string("H+");
    
    // Info.
    //
    info_subphysics_ = std::string("ICDOWSubPhysics");
}

ICDOWSubPhysics::~ICDOWSubPhysics(){
//    delete wavecurvefactory_;

//    delete odesolver_;

//    delete compositecurve_;

//    delete shockcurve_;

//    for (int i = 0; i < extension_curve.size(); i++) delete extension_curve[i];

    delete inflection_curve_;

    delete rarefactioncurve_;

    delete coincidence_contour_;
    delete coincidence_;

//    // Not sure if this should really be done like this.
//    // Perhaps it is best to eliminate only the HugoniotCurves that were instantiated
//    // in this class, and let the rest be deleted at the father's dtor. 
//    //
//    for (int i = 0; i < hugoniot_curve.size(); i++) delete hugoniot_curve[i];

    delete gridvalues_;

    delete boundary_;
    delete accumulation_;
    delete flux_;

    delete chemistry_;
    delete hydrodynamics_;

    for (int i = equation_parameter_.size() - 1; i >= 0; i--) delete equation_parameter_[i];
}


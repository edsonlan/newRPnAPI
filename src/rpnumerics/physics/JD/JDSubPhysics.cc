#include "JDSubPhysics.h"

JDSubPhysics::JDSubPhysics() : SubPhysics() {
    // Flux, accumulation.
    //
    Parameter *epsilon_parameter = new Parameter(std::string("epsilon"), 1e-1);
    equation_parameter_.push_back(epsilon_parameter);

    flux_         = new JDFluxFunction(epsilon_parameter);
    accumulation_ = new JDAccumulationFunction;
    
    // Boundary.
    RealVector pmin(2), pmax(2);
    pmin(0) = 0.1; // THERE IS A PROBLEM IF u = 0!!!
    pmin(1) = 0.0;

    pmax = pmin + 5.0;
    boundary_ = new RectBoundary(pmin, pmax);
    
    // GridValues.
    //
    std::vector<int> number_of_cells(2);
    number_of_cells[0] = 256;
    number_of_cells[1] = 256;

    gridvalues_ = new GridValues(boundary_, pmin, pmax, number_of_cells);

    // Implicit Hugoniot.
    //
    ImplicitHugoniotCurve *ihc = new ImplicitHugoniotCurve(flux_, accumulation_, boundary_);
    ihc->subphysics(this);

    hugoniot_curve.push_back(ihc);

    // Coincidence.
    //
    coincidence_ = new CoincidenceJD((const JDFluxFunction*)flux_, (const JDAccumulationFunction*)accumulation_);
    coincidence_contour_ = new Coincidence_Contour(coincidence_);

    // Evaporation extension.
    //
    evapext = new JDEvap_Extension((const JDFluxFunction*)flux_, cjd);

    // RarefactionCurve.
    //
    rarefactioncurve_ = new RarefactionCurve(accumulation_, flux_, boundary_);

    // HugoniotContinuation.
    //
    hugoniotcontinuation_ = new HugoniotContinuation2D2D(flux_, accumulation_, boundary_);
    
    hugoniot_curve.push_back(hugoniotcontinuation_);

    // ShockCurve.
    //
    shockcurve_ = new ShockCurve(hugoniotcontinuation_);

    // Composite.
    //
    compositecurve_ = new JDEvaporationCompositeCurve(accumulation_, flux_, boundary_, evapext);

    // ODE_Solver.
    //
    odesolver_ = new LSODE;

    // WaveCurveFactory.
    //
    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);

    // Inflection.
    //
    inflection_curve_ = new Inflection_Curve;

    // Canvas-related.
    //
    transformation_matrix_ = DoubleMatrix::eye(2);

    xlabel_ = std::string("u");
    ylabel_ = std::string("v");

    info_subphysics_ = std::string("JDSubPhysics");
}

JDSubPhysics::~JDSubPhysics(){
    delete inflection_curve_;

    delete wavecurvefactory_;
    delete odesolver_;
    delete compositecurve_;
    delete shockcurve_;
//    delete hugoniotcontinuation_;
    delete rarefactioncurve_;

    delete evapext;
    delete coincidence_contour_;
    delete coincidence_;

    // Not sure if this should really be done like this.
    // Perhaps it is best to eliminate only the HugoniotCurves that were instantiated
    // in this class, and let the rest be deleted at the father's dtor. 
    //
    for (int i = 0; i < hugoniot_curve.size(); i++) delete hugoniot_curve[i];
    
    delete gridvalues_;
    delete boundary_;
    delete accumulation_;
    delete flux_;
    
    for (int i = equation_parameter_.size() - 1; i >=0; i--) delete equation_parameter_[i];
    

}

void JDSubPhysics::shock_cases(std::vector<int> &type, std::vector<std::string> &name) const {
    type.clear();
    name.clear();

    type.push_back(JD_GENERIC_POINT);
    name.push_back(std::string("Generic point"));

    return;
}

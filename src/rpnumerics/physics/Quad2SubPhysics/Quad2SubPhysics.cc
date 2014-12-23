#include "Quad2SubPhysics.h"

Quad2SubPhysics::Quad2SubPhysics(): SubPhysics(){
    // Parameters.
    //
    double a = -1.0;
    double b =  0.0;
    double c =  0.1;

    a1_parameter = new Parameter(std::string("a1"), a);
    b1_parameter = new Parameter(std::string("b1"), b);
    c1_parameter = new Parameter(std::string("c1"), 1.0);
    d1_parameter = new Parameter(std::string("d1"), 0.0);
    e1_parameter = new Parameter(std::string("e1"), c);

    a2_parameter = new Parameter(std::string("a2"), b);
    b2_parameter = new Parameter(std::string("b2"), 1.0);
    c2_parameter = new Parameter(std::string("c2"), 0.0);
    d2_parameter = new Parameter(std::string("d2"), -c);
    e2_parameter = new Parameter(std::string("e2"), 0.0);

    equation_parameter_.push_back(a1_parameter);
    equation_parameter_.push_back(b1_parameter);
    equation_parameter_.push_back(c1_parameter);
    equation_parameter_.push_back(d1_parameter);
    equation_parameter_.push_back(e1_parameter);
    equation_parameter_.push_back(a2_parameter);
    equation_parameter_.push_back(b2_parameter);
    equation_parameter_.push_back(c2_parameter);
    equation_parameter_.push_back(d2_parameter);
    equation_parameter_.push_back(e2_parameter);

    flux_ = new Quad2FluxFunction(a1_parameter, b1_parameter, c1_parameter, d1_parameter, e1_parameter,
                                  a2_parameter, b2_parameter, c2_parameter, d2_parameter, e2_parameter);

    accumulation_ = new Quad2AccumulationFunction;

    // Boundary.
    //
    RealVector pmin(2), pmax(2);
    pmin(0) = pmin(1) = -.5;
    pmax = -pmin;

    boundary_ = new RectBoundary(pmin, pmax);

    // GridValues.
    //
    std::vector<int> number_of_cells(2);
    number_of_cells[0] = 128;
    number_of_cells[1] = 128;

    gridvalues_ = new GridValues(boundary_, boundary_->minimums(), boundary_->maximums(), number_of_cells);

    // Rarefaction.
    //
    rarefactioncurve_ = new RarefactionCurve(accumulation_, flux_, boundary_);

    // Explicit Hugoniot.
    //
    Quad2ExplicitHugoniotCurve *ehc = new Quad2ExplicitHugoniotCurve(this);
    hugoniot_curve.push_back(ehc);

    // Implicit Hugoniot.
    //
    ImplicitHugoniotCurve *ihc = new ImplicitHugoniotCurve(flux_, accumulation_, boundary_);
    ihc->subphysics(this);

    hugoniot_curve.push_back(ihc);

    // HugoniotContinuation.
    //
    hugoniotcontinuation_ = new HugoniotContinuation2D2D(flux_, accumulation_, boundary_);

    hugoniot_curve.push_back(hugoniotcontinuation_);

    // Shockcurve.
    //
    shockcurve_ = new ShockCurve(hugoniotcontinuation_);

    // Composite.
    //
    compositecurve_ = new CompositeCurve(accumulation_, flux_, boundary_, shockcurve_, 0/*&bc*/);

    odesolver_ = new LSODE;

    // Wavecurvefactory.
    //
    wavecurvefactory_ = new WaveCurveFactory(accumulation_, flux_, boundary_, odesolver_, rarefactioncurve_, shockcurve_, compositecurve_);

    // Inflection.
    //
    inflection_curve_ = new Inflection_Curve;

    // Double_Contact.
    //
    doublecontact_ = new Double_Contact;

    // Extension curve.
    //
    Implicit_Extension_Curve *iec = new Implicit_Extension_Curve;
    extension_curve.push_back(iec);

    // Canvas-related.
    //
    transformation_matrix_ = DoubleMatrix::eye(2);

    xlabel_ = std::string("u");
    ylabel_ = std::string("v");

    info_subphysics_ = std::string("Quad2SubPhysics");
}

Quad2SubPhysics::~Quad2SubPhysics(){
    for (int i = extension_curve.size() - 1; i >= 0; i--) delete extension_curve[i];

    delete doublecontact_;
    delete inflection_curve_;
    delete wavecurvefactory_;
    delete odesolver_;
    delete compositecurve_;
    delete shockcurve_;

    for (int i = hugoniot_curve.size() - 1; i >= 0; i--) delete hugoniot_curve[i];

    delete rarefactioncurve_;
    delete gridvalues_;
    delete boundary_;
    delete accumulation_;
    delete flux_;

    for (int i = 0; i < equation_parameter_.size(); i++) delete equation_parameter_[i];
}


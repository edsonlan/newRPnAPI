#include "TPCWSubPhysics.h"

TPCWSubPhysics::TPCWSubPhysics(){
    // Flux parameters.
    //
    abs_perm_parameter = new Parameter(std::string("abs_perm"), 3e-12);
    sin_beta_parameter = new Parameter(std::string("sin_beta"), 0.0);
    cnw_parameter      = new Parameter(std::string("cnw"), 0.0);
    cng_parameter      = new Parameter(std::string("cng"), 0.0);
    expw_parameter     = new Parameter(std::string("expw"), 2.0);
    expg_parameter     = new Parameter(std::string("expg"), 2.0);

    equation_parameter_.push_back(abs_perm_parameter);
    equation_parameter_.push_back(sin_beta_parameter);
    equation_parameter_.push_back(cnw_parameter);
    equation_parameter_.push_back(cng_parameter);
    equation_parameter_.push_back(expw_parameter);
    equation_parameter_.push_back(expg_parameter);

    // Accumulation parameter.
    //
    phi_parameter = new Parameter(std::string("phi"), 0.38);

    equation_parameter_.push_back(phi_parameter);

    // Molar density.
    //
    mdv = new MolarDensity(MOLAR_DENSITY_VAPOR,  100.900000e5);
    mdl = new MolarDensity(MOLAR_DENSITY_LIQUID, 100.900000e5);

    flash = new VLE_Flash_TPCW(mdl, mdv);

    tc = new Thermodynamics("../../../../../c++/rpnumerics/physics/CompositionalPhysics/TPCWPhysics/TPCWSubPhysics/hsigmaC_spline.txt");
    tc->set_flash(flash);

    bool has_gravity = false;
    bool has_horizontal = true;

    flux_ = new Flux2Comp2PhasesAdimensionalized(abs_perm_parameter, sin_beta_parameter, 
                                                 cnw_parameter, cng_parameter,
                                                 expw_parameter, expg_parameter,
                                                 has_gravity, has_horizontal,
                                                 tc);

    accumulation_ = new Accum2Comp2PhasesAdimensionalized(phi_parameter, tc);

    // Boundary.
    //
    double Theta_min = 0.099309;
    double Theta_max = 0.576511;

    RealVector pmin(3), pmax(3);

    pmin.component(0) = 0.0;
    pmin.component(1) = Theta_min;
    pmin.component(2) = 0.0;

    pmax.component(0) = 1.0;
    pmax.component(1) = Theta_max;
    pmax.component(2) = 2.0;

    boundary_     = new RectBoundary(pmin, pmax);

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

    // Canvas-related.
    //
    transformation_matrix_ = DoubleMatrix::eye(3);
    transformation_matrix_(0, 0) = -1.0;

    xlabel_ = std::string("s");
    ylabel_ = std::string("Theta");
    
    // Info.
    //
    info_subphysics_ = std::string("TPCWSubPhysics");
}

TPCWSubPhysics::~TPCWSubPhysics(){
    // Not sure if this should really be done like this.
    // Perhaps it is best to eliminate only the HugoniotCurves that were instantiated
    // in this class, and let the rest be deleted at the father's dtor. 
    //
    for (int i = 0; i < hugoniot_curve.size(); i++) delete hugoniot_curve[i];

    delete gridvalues_;
    delete boundary_;
    delete accumulation_;
    delete flux_;
    delete tc;
    delete flash;
    delete mdl;
    delete mdv;

    for (int i = equation_parameter_.size() - 1; i >= 0; i--) delete equation_parameter_[i];
}


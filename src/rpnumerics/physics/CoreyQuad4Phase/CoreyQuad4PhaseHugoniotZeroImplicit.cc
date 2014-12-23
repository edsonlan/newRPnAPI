#include "CoreyQuad4PhaseHugoniotZeroImplicit.h"

RealVector CoreyQuad4PhaseHugoniotZeroImplicit::sw_zero(const RealVector &p){
    RealVector pp(3);

    pp(0) = 0.0;
    pp(1) = p(0);
    pp(2) = p(1);

    return pp;
}

RealVector CoreyQuad4PhaseHugoniotZeroImplicit::so_zero(const RealVector &p){
    RealVector pp(3);

    pp(0) = p(0);
    pp(1) = 0.0;
    pp(2) = p(1);

    return pp;
}

RealVector CoreyQuad4PhaseHugoniotZeroImplicit::sg_zero(const RealVector &p){
    RealVector pp(3);

    pp(0) = p(0);
    pp(1) = p(1);
    pp(2) = 0.0;
}

RealVector CoreyQuad4PhaseHugoniotZeroImplicit::sc_zero(const RealVector &p){
    RealVector pp(3);

    pp(0) = p(0);
    pp(1) = p(1);
    pp(2) = 1.0 - p(0) - p(1);
}

RealVector CoreyQuad4PhaseHugoniotZeroImplicit::f_zero(const RealVector &p){
    RealVector pp = (*convert)(p);

    JetMatrix Fjet, Gjet;
    flux->jet(pp, Fjet, 0);
    // accumulation->jet(pp, Gjet, 0);

    RealVector f(2);

    RealVector dflux  = Fjet.function() - ref.F;
    RealVector daccum = pp - ref.point;

    f(0) = daccum(1)*dflux(0) - daccum(0)*dflux(1);
    f(1) = daccum(2)*dflux(1) - daccum(1)*dflux(2);

    return f;
}

CoreyQuad4PhaseHugoniotZeroImplicit::CoreyQuad4PhaseHugoniotZeroImplicit(): ZeroImplicit(){
}

CoreyQuad4PhaseHugoniotZeroImplicit::~CoreyQuad4PhaseHugoniotZeroImplicit(){
}

void CoreyQuad4PhaseHugoniotZeroImplicit::find_zeroes(int side, std::vector<RealVector> &zeroes){
    if      (side == SW_ZERO) convert = &sw_zero;
    else if (side == SO_ZERO) convert = &so_zero;
    else if (side == SG_ZERO) convert = &sg_zero;
    else if (side == SC_ZERO) convert = &sc_zero;

    std::vector<RealVector> temp_points;

    // Map back to 3D.
    //
    zeroes.clear();
    for (int i = 0; i < temp_points.size(); i++) zeroes.push_back((*convert)(temp_points[i]));

    return;
}


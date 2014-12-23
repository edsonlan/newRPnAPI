#include "CoreyQuad_Params.h"

CoreyQuad_Params::CoreyQuad_Params(const double grw, const double grg, const double gro,
                           const double muw, const double mug, const double muo,
                           const double vel
                           ) : FluxParams(RealVector(7)){
    component(0, grw);
    component(1, grg);
    component(2, gro);

    component(3, muw);
    component(4, mug);
    component(5, muo);

    component(6, vel);

   
}

CoreyQuad_Params::CoreyQuad_Params() : FluxParams(RealVector(7)){
    component(0, 0.0);
    component(1, 0.0);
    component(2, 0.0);

    component(3, 1.0);
    component(4, 1.0);
    component(5, 1.0);

    component(6, 1.0);

  
}

CoreyQuad_Params::CoreyQuad_Params(const RealVector & params):FluxParams(params){}
CoreyQuad_Params::CoreyQuad_Params(const CoreyQuad_Params &copy):FluxParams(copy.params()){}




CoreyQuad_Params::~CoreyQuad_Params(){
}

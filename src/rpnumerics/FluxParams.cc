
#include "FluxParams.h"


FluxParams::FluxParams(const RealVector & params) {
    
    params_= new RealVector(params.size());
    
    for (int i=0; i< params_->size();i++){
        
        params_->component(i)=params.component(i);
    }
    
    
}

FluxParams::FluxParams(const int size,  double * coords):params_(new RealVector(size, coords)){}

FluxParams::FluxParams(const FluxParams &params){
    
    params_= new RealVector(params.params());
    
    int i ;
    
    for (i=0;i < params_->size();i++){
        
        params_->component(i)=params.component(i);
        
    }
}


FluxParams::~FluxParams() {
    delete params_;
    
}



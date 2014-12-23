#include "AccumulationFunction.h"

AccumulationFunction::AccumulationFunction(void) : params_(new AccumulationParams()) { }



//RpFunction * AccumulationFunction::clone() const  { return new AccumulationFunction(*this);}


AccumulationFunction::AccumulationFunction(const AccumulationFunction & copy ){
     params_= new AccumulationParams(copy.accumulationParams().params());
 }


//int AccumulationFunction::jet(const WaveState &u, JetMatrix &m, int degree) const{
//
//    return 0;
//}


AccumulationFunction::~AccumulationFunction(void){delete params_;}


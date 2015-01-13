#include "RPnMethod.h"

RPnMethod::RPnMethod(Configuration * configuration) : config_(configuration) {
}

const Configuration & RPnMethod::getConfiguration() const {
    return *config_;
}




RPnMethod::~RPnMethod() {
   
    delete config_;
}
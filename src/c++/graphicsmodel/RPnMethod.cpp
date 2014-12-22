#include "RPnMethod.h"

RPnMethod::RPnMethod( Configuration * configuration) : config_(configuration) {
}

const Configuration & RPnMethod::getConfiguration() const {
    return *config_;
}

 Configuration & RPnMethod::getConfiguration()  {
    return *config_;
}


RPnMethod::~RPnMethod() {
    cout<<"Chamando destrutor do metodo"<<endl;
    delete config_;
}
#include "Quad4AccumulationFunction.h"

Quad4AccumulationFunction::Quad4AccumulationFunction(void) {
}

Quad4AccumulationFunction * Quad4AccumulationFunction::clone() const {
    return new Quad4AccumulationFunction(*this);
}



Quad4AccumulationFunction::~Quad4AccumulationFunction(void) {}

int Quad4AccumulationFunction::jet(const WaveState&, JetMatrix&, int) const {
    return 0;
}










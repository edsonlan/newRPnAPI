#include "Quad3AccumulationFunction.h"
#include <stdio.h>

Quad3AccumulationFunction::Quad3AccumulationFunction(void) {
}

Quad3AccumulationFunction * Quad3AccumulationFunction::clone() const {
    return new Quad3AccumulationFunction(*this);
}



Quad3AccumulationFunction::~Quad3AccumulationFunction(void) {}

int Quad3AccumulationFunction::jet(const WaveState&, JetMatrix&, int) const {
    return 0;
}










#include "QuadWaveState.h"


QuadWaveState::QuadWaveState(const int dim) : 
	WaveState(dim){}

QuadWaveState::QuadWaveState(const int dim, sigma sig) : 
	WaveState(dim),sigma_(sig){}

QuadWaveState::QuadWaveState(const RealVector & vector) :
	WaveState(vector){}

QuadWaveState::QuadWaveState(const RealVector & vector, const sigma sig) :
	WaveState(vector),sigma_(sig){}

QuadWaveState::~QuadWaveState(){}



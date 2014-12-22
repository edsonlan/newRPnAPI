/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) WaveState.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "WaveState.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */

WaveState::WaveState(const RealVector & coords):coords_(new RealVector(coords)){}

WaveState::WaveState(const int dim ):coords_(new RealVector(dim)) {}


WaveState::~WaveState() {
    delete coords_;
}

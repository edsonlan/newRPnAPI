/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) ShockFlowParams.cc
 */
#include "ShockFlowParams.h"

ShockFlowParams::ShockFlowParams(const PhasePoint &p, double sigma) :  phasePoint_(new PhasePoint(p)),sigma_(sigma) {
}

ShockFlowParams::ShockFlowParams(const PhasePoint &p) : phasePoint_(new PhasePoint(p)),sigma_(0.0) {
}

ShockFlowParams::ShockFlowParams(const ShockFlowParams & copy) :  phasePoint_(new PhasePoint(copy.getPhasePoint())),sigma_(copy.getSigma()){}




ShockFlowParams::~ShockFlowParams() {
    delete phasePoint_;
}

const PhasePoint & ShockFlowParams::getPhasePoint() const {

    return *phasePoint_;
}

void ShockFlowParams::setPhasePoint(const PhasePoint &p) {

    delete phasePoint_;
    phasePoint_ = new PhasePoint(p);

}


double ShockFlowParams::getSigma() const {

    return sigma_;
}

void ShockFlowParams::setSigma(double sigma) {

    sigma_ = sigma;
}


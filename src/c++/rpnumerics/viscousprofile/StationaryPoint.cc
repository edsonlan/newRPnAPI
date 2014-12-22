/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) StationaryPoint.cc
 */

/*
 * ---------------------------------------------------------------
 * Includes:
 */
#include "StationaryPoint.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */


StationaryPoint::StationaryPoint(const RealVector & coords, int type) : coords_(new RealVector(coords)),type_(type) {
}

const RealVector & StationaryPoint::coords(){return *coords_;}
int StationaryPoint::type(){return type_;}



StationaryPoint::~StationaryPoint() {
    delete coords_;

}
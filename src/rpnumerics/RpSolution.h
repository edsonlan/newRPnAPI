/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RpSolution.h
 */

#ifndef _RpSolution_H
#define _RpSolution_H

/*
 * ---------------------------------------------------------------
 * Includes:
 */


/*
 * ---------------------------------------------------------------
 * Definitions:
 */


class RpSolution {
    
private:
    
public:
    static const int DEFAULT_NULL_FLAG = -1;
    static const int OUT_OF_THE_BOUNDARY = 0;
    static const int INTERSECTION_WITH_POINCARE = 1;
    static const int MAX_NUM_OF_ITERATION_REACHED = 2;
    // Stationary point calculations
    static const int NO_CONVERGENCE_IN_STATIONARY_POINT_COMPUTATION = 3;
    static const int CONVERGENCE_IN_STATIONARY_POINT_COMPUTATION = 4;
    static const int FOUND_OUT_OF_BOUNDARY = 5;
    // Manifold calculations
    static const int MANIFOLD_DOES_NOT_EXIST = 6;
    static const int MANIFOLD_FINISHES_ON_POINCARE_SECTION = 7;
    static const int MANIFOLD_FINISHES_ON_BOUNDARY = 8;
    static const int MANIFOLD_FINISHES_DUE_TO_MANY_STEPS = 9;
    static const int STATIONARY_POINT_IS_TOO_CLOSE_TO_POINCARE_SECTION = 10;
    static const int STATIONARY_POINT_IS_TOO_CLOSE_TO_BOUNDARY = 11;
    // Connections calculations
    static const int CONNECTION_NOT_FINISHED = 12;
    static const int DEVIATION_INCREASED = 13;
    static const int NO_POINCARE_SECTION = 14;
    static const int CONNECTED = 15;
    
};

#endif //! _RpSolution_H


































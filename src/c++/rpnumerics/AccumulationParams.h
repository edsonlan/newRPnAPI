/* IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) AccumulationParams.h
 **/

#ifndef _AccumulationParams_H
#define	_AccumulationParams_H

#include "RealVector.h"

class AccumulationParams {
private:

    RealVector * params_;


public:

    AccumulationParams();
    AccumulationParams(const RealVector & params);
    AccumulationParams(const int size, double *);
    AccumulationParams(const AccumulationParams &);
    virtual ~AccumulationParams();

    const RealVector & params(void) const;
    const double component(int index) const;
    void params(int size, double * coords);
    void params(const RealVector & params);
    void component(int index, double value);

    AccumulationParams defaultParams(void);

};

inline const RealVector & AccumulationParams::params(void) const {
    return *params_;
}

inline const double AccumulationParams::component(int index) const {
    return params_->component(index);
}

inline void AccumulationParams::params(int size, double * coords) {
    delete params_;
    params_ = new RealVector(size, coords);

}

inline void AccumulationParams::params(const RealVector & params) {
    delete params_;
    params_ = new RealVector(params);
}

inline void AccumulationParams::component(int index, double value) {
    params_->component(index) = value;
}



#endif	/* _AccumulationParams_H */


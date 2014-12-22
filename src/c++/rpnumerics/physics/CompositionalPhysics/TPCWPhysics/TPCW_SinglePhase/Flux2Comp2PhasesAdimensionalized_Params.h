#ifndef _FLUX2COMP2PHASESADIMENSIONALIZED_PARAMS_
#define _FLUX2COMP2PHASESADIMENSIONALIZED_PARAMS_

#include "FluxParams.h"
#include <iostream>

#include "Thermodynamics.h"


using namespace std;

class Flux2Comp2PhasesAdimensionalized_Params : public FluxParams {
private:
    double abs_perm; // = 3e-12
    double sin_beta; // = 0.0
    const double const_gravity; //= 9.8;
    bool has_gravity_, has_horizontal_;

    Thermodynamics *TD_;

protected:
public:


    Flux2Comp2PhasesAdimensionalized_Params(const RealVector &,Thermodynamics *);

    Flux2Comp2PhasesAdimensionalized_Params(double abs_perm, double sin_beta, double const_gravity,
            bool has_gravity,
            bool has_horizontal,
            Thermodynamics * TD);


    Flux2Comp2PhasesAdimensionalized_Params(const Flux2Comp2PhasesAdimensionalized_Params &);



    virtual ~Flux2Comp2PhasesAdimensionalized_Params();

    Thermodynamics * get_thermodynamics(void) const;

    bool has_gravity(void) const;
    bool has_horizontal(void) const;
};

#endif // _FLUX2COMP2PHASESADIMENSIONALIZED_PARAMS_


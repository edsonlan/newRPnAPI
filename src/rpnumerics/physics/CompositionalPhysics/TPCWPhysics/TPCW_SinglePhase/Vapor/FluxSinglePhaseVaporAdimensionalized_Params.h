#ifndef _FLUXSINGLEPHASEVAPORADIMENSIONALIZED_PARAMS_
#define _FLUXSINGLEPHASEVAPORADIMENSIONALIZED_PARAMS_

#include <string>

#include "FluxParams.h"
//#include "FluxSinglePhaseVaporAdimensionalized.h"

//#include "Thermodynamics_Common.h"

class FluxSinglePhaseVaporAdimensionalized_Params : public FluxParams {
    private:
        double Tref_rock, Tref_water, P;
        double Rock_Cr;
        double Cw;
        double T_typical_;
        double Rho_typical_;
        double U_typical_;
        double h_typical_;

        std::string hsigmaC;
    protected:
    public:
        FluxSinglePhaseVaporAdimensionalized_Params(double mc, double mw, double Press, 
                                                    const char *hsigmaC,
                                                    double Tr, double Tw,
                                                    double Cr,
                                                    double cw,
                                                    double T_typical,
                                                    double Rho_typical,
                                                    double U_typical,
                                                    double h_typical);
        ~FluxSinglePhaseVaporAdimensionalized_Params();

        const char * get_name(void);
};

#endif // _FLUXSINGLEPHASEVAPORADIMENSIONALIZED_PARAMS_


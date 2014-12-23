#ifndef _FLUX2COMP2PHASESADIMENSIONALIZED_
#define _FLUX2COMP2PHASESADIMENSIONALIZED_

#include <stdio.h>
#include <stdlib.h>
#include "FluxFunction.h"

#include "Thermodynamics.h"
#include "Flux2Comp2PhasesAdimensionalized_Params.h"

#define FLUX2COMP2PHASESADIMENSIONALIZED_PURE_GRAVITY 0
#define FLUX2COMP2PHASESADIMENSIONALIZED_GRAVITY      1
#define FLUX2COMP2PHASESADIMENSIONALIZED_HORIZONTAL   2

#define JETTESTER_ENABLED_TPCWFLUX


class Flux2Comp2PhasesAdimensionalized: public FluxFunction
                                       
{
    private:
    protected:
        class FracFlow2PhasesVerticalAdimensionalized {
            private:
                Flux2Comp2PhasesAdimensionalized * fluxComplete_;
            public:
                FracFlow2PhasesVerticalAdimensionalized(Flux2Comp2PhasesAdimensionalized *);
                virtual ~FracFlow2PhasesVerticalAdimensionalized(){}

                int Diff_FracFlow2PhasesVerticalAdimensionalized(double sw, double Theta, int degree, JetMatrix &m);
        };

        class FracFlow2PhasesHorizontalAdimensionalized {
            private:
                Flux2Comp2PhasesAdimensionalized * fluxComplete_;
            public:
                FracFlow2PhasesHorizontalAdimensionalized(Flux2Comp2PhasesAdimensionalized *);
                virtual ~FracFlow2PhasesHorizontalAdimensionalized(){}

                int Diff_FracFlow2PhasesHorizontalAdimensionalized(double sw, double Theta, int degree, JetMatrix &m);
        };

        class ReducedFlux2Comp2PhasesAdimensionalized {
            private:
                Flux2Comp2PhasesAdimensionalized * fluxComplete_;
            public:
                ReducedFlux2Comp2PhasesAdimensionalized(Flux2Comp2PhasesAdimensionalized *);
                virtual ~ReducedFlux2Comp2PhasesAdimensionalized(){}

                int jet(const WaveState &u, JetMatrix &m, int degree) const;
        };

        // Fluid dynamics
        bool has_gravity;
        bool has_horizontal;
        double const_gravity;

        Parameter *abs_perm_parameter, *sin_beta_parameter;
        Parameter *cnw_parameter,  *cng_parameter;
        Parameter *expw_parameter, *expg_parameter;

        // Thermodynamics
        Thermodynamics *TD;

        // FracFlows
        FracFlow2PhasesHorizontalAdimensionalized *FH;
        FracFlow2PhasesVerticalAdimensionalized *FV;

        //ReducedFlux
        ReducedFlux2Comp2PhasesAdimensionalized *reducedFlux;

        double T_typical_;
    public:
        Flux2Comp2PhasesAdimensionalized(Parameter *abs_perm, Parameter *sin_beta, 
                                         Parameter *cnw, Parameter *cng,
                                         Parameter *expw, Parameter *expg,
                                         bool has_grav, bool has_hor,
                                         Thermodynamics *td);
        virtual ~Flux2Comp2PhasesAdimensionalized();

        Thermodynamics* getThermo()const ;

        FracFlow2PhasesHorizontalAdimensionalized * getHorizontalFlux()const ;
        FracFlow2PhasesVerticalAdimensionalized * getVerticalFlux()const;
        ReducedFlux2Comp2PhasesAdimensionalized * getReducedFlux()const;

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
        void type(int t);

        const FracFlow2PhasesHorizontalAdimensionalized* get_horizontal_flux(){return FH;}
        const FracFlow2PhasesVerticalAdimensionalized*   get_vertical_flux(){return FV;}

        #ifdef JETTESTER_ENABLED_TPCWFLUX
        static int testable_flux_jet(void *obj, const RealVector &state, int degree, JetMatrix &jm){
            const Flux2Comp2PhasesAdimensionalized *flux = (const Flux2Comp2PhasesAdimensionalized*)obj;

            int info = flux->jet(state, jm, degree);
            return info;
        }

        void list_of_functions(std::vector<int (*)(void*, const RealVector&, int degree, JetMatrix&)> &list,
                               std::vector<std::string> &name){
            list.clear();
            list.push_back(&testable_flux_jet);

            name.clear();
            name.push_back(std::string("TPCW Flux"));

            return;
        }
        #endif
};

#endif // _FLUX2COMP2PHASESADIMENSIONALIZED_


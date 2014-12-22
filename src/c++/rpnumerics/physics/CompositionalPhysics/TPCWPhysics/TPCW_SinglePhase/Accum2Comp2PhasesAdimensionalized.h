#ifndef _ACCUM2COMP2PHASESADIMENSIONALIZED_
#define _ACCUM2COMP2PHASESADIMENSIONALIZED_

#include <stdio.h>
#include <stdlib.h>
#include "AccumulationFunction.h"

#include "Thermodynamics.h"
#include "Parameter.h"

#define JETTESTER_ENABLED_TPCWACCUMULATION


class Accum2Comp2PhasesAdimensionalized: public AccumulationFunction
{
    private:
    protected:
        class ReducedAccum2Comp2PhasesAdimensionalized:public AccumulationFunction {
        private:
            Accum2Comp2PhasesAdimensionalized * totalAccum_;
        public:
            ReducedAccum2Comp2PhasesAdimensionalized(Accum2Comp2PhasesAdimensionalized *);
            virtual ~ReducedAccum2Comp2PhasesAdimensionalized();

            int jet(const WaveState &u, JetMatrix &m, int degree) const; // ONLY DEGREE 0
        };

        Parameter *phi_parameter_;
        Thermodynamics *TD;
        ReducedAccum2Comp2PhasesAdimensionalized * reducedAccum_;
    public:
        Accum2Comp2PhasesAdimensionalized(Parameter *phi, Thermodynamics *td);

        virtual ~Accum2Comp2PhasesAdimensionalized();

        int jet(const WaveState &u, JetMatrix &m, int degree) const;
        ReducedAccum2Comp2PhasesAdimensionalized* getReducedAccumulation() const;

        #ifdef JETTESTER_ENABLED_TPCWACCUMULATION
        static int testable_accumulation_jet(void *obj, const RealVector &state, int degree, JetMatrix &jm){
            const Accum2Comp2PhasesAdimensionalized *accum = (const Accum2Comp2PhasesAdimensionalized*)obj;

            int info = accum->jet(state, jm, degree);
            return info;
        }

        void list_of_functions(std::vector<int (*)(void*, const RealVector&, int degree, JetMatrix&)> &list,
                               std::vector<std::string> &name){
            list.clear();
            list.push_back(&testable_accumulation_jet);

            name.clear();
            name.push_back(std::string("TPCW Accumulation"));

            return;
        }
        #endif  
};

#endif // _ACCUM2COMP2PHASESADIMENSIONALIZED_


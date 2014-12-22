#ifndef _ICDOWACCUMULATIONFUNCTION_
#define _ICDOWACCUMULATIONFUNCTION_

#include "AccumulationFunction.h"
#include "ICDOWChemistry.h"

class ICDOWAccumulationFunction: public AccumulationFunction
                                 
{
    private:
    protected:
        ICDOWChemistry *chemistry;
        
        Parameter *phi_parameter;
    public:
        ICDOWAccumulationFunction(Parameter *phi, ICDOWChemistry *ch);
        virtual ~ICDOWAccumulationFunction();

        int reduced_jet(const WaveState &state, JetMatrix &m, int degree) const;
        int reduced_jet(const RealVector &u, JetMatrix &m, int degree) const {
            WaveState w(u);

            int info = reduced_jet(w, m, degree);

            return info;
        }

        int jet(const WaveState &u, JetMatrix &m, int degree) const;

       
};

#endif // _ICDOWACCUMULATIONFUNCTION_


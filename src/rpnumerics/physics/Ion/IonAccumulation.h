#ifndef _IONACCUMULATION_
#define _IONACCUMULATION_

#include "AccumulationFunction.h"
#include "IonAdsorption.h"

class IonAccumulation : public AccumulationFunction {
    private:
    protected:
        const IonAdsorption *adsorption;
    public:
        IonAccumulation(const IonAdsorption *a);
        virtual ~IonAccumulation();

        int jet(const WaveState &w, JetMatrix &m, int degree) const;
        IonAccumulation * clone() const;
};

#endif // _IONACCUMULATION_


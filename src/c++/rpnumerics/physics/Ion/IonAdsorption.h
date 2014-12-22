#ifndef _IONADSORPTION_
#define _IONADSORPTION_

// TODO: This class will be derived from AuxiliaryClass in the future!

#include "WaveState.h"
#include "JetMatrix.h"

class IonAdsorption {
    private:
    protected:
    public:
        IonAdsorption();
        virtual ~IonAdsorption();

        int jet(const WaveState &w, JetMatrix &a, int degree) const;
};

#endif // _IONADSORPTION_


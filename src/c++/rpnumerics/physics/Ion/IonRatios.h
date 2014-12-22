#ifndef _IONRATIOS_
#define _IONRATIOS_

// TODO: This class will be derived from AuxiliaryClass in the future!

#include "WaveState.h"
#include "JetMatrix.h"

class IonRatios {
    private:
    protected:
    public:
        IonRatios();
        virtual ~IonRatios();

        int jet(const WaveState &w, JetMatrix &r, int degree) const;
};

#endif // _IONRATIOS_


#ifndef _IONPERMEABILITY_
#define _IONPERMEABILITY_

// TODO: This class will be derived from AuxiliaryClass in the future!

#include "WaveState.h"
#include "JetMatrix.h"

class IonPermeability {
    private:
    protected:
    public:
        IonPermeability();
        virtual ~IonPermeability();

        int jet(const WaveState &w, JetMatrix &p, int degree) const;
};

#endif // _IONPERMEABILITY_


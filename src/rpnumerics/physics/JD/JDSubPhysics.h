#ifndef _JDSUBPHYSICS_
#define _JDSUBPHYSICS_

#include "SubPhysics.h"
#include "JDFluxFunction.h"
#include "JDAccumulationFunction.h"
#include "RectBoundary.h"
#include "CoincidenceJD.h"
#include "JDEvap_Extension.h"
#include "JDEvaporationCompositeCurve.h"
#include "HugoniotContinuation2D2D.h"
#include "LSODE.h"

#define JD_GENERIC_POINT 0

class JDSubPhysics : public SubPhysics {
    private:
    protected:
        Parameter *epsilon_parameter;
            
        CoincidenceJD *cjd;
        JDEvap_Extension *evapext;
    public:
        JDSubPhysics();
        virtual ~JDSubPhysics();
        
        void shock_cases(std::vector<int> &type, std::vector<std::string> &name) const;
};

#endif // _JDSUBPHYSICS_

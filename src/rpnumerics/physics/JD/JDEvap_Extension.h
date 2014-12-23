#ifndef _JDEVAP_EXTENSION_
#define _JDEVAP_EXTENSION_

#include "Extension.h"

#include "Coincidence.h"
#include "JDFluxFunction.h"
#include "Utilities.h"

// Only valid for TPCW-like classes.
//
class JDEvap_Extension : public Extension {
    private:
    protected:
        // The component index such that p(i) = extension(p(i)).
        //
        unsigned int index_of_constant;

        // The other component.
        //
        unsigned int index_of_non_constant;

        const Coincidence *coincidence_;

        const JDFluxFunction *flux_;
    public:
        JDEvap_Extension(const JDFluxFunction *f, const Coincidence *c);
        virtual ~JDEvap_Extension();

        // 
        int extension(const RealVector &p, RealVector &ext_p);

        std::string name() const {
            return std::string("JD Evaporation extension");
        }

        int extension_type(){return EXPLICIT_EXTENSION;}

        // This method is only valid for this extension (so far).
        const Coincidence * coincidence() const {return coincidence_;}
};

#endif // _JDEVAP_EXTENSION_


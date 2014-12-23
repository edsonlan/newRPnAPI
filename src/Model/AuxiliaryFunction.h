#ifndef _AUXILIARYFUNCTION_
#define _AUXILIARYFUNCTION_

#include <vector>
#include "Parameter.h"
#include "JetMatrix.h"

class AuxiliaryFunction {
    private:
        // Users should not add any private members.
    protected:
        std::string info_auxiliaryfunction_;

        std::vector<Parameter*> parameters_;
    public:
        AuxiliaryFunction(){}
        virtual ~AuxiliaryFunction(){}

        virtual std::string info_auxiliary_function(){return info_auxiliaryfunction_;}
        virtual void parameter(std::vector<Parameter*> &vp){vp = parameters_; return;}
};

#endif // _AUXILIARYFUNCTION_


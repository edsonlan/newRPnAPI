#ifndef _THREEIMPLICITFUNCTIONS_
#define _THREEIMPLICITFUNCTIONS_

#include <vector>
#include "GridValues.h"

class Contour2x2_Method;

class ThreeImplicitFunctions {
    private:
    protected:
        GridValues *gv_left, *gv_right;
        bool singular;
    public:
        ThreeImplicitFunctions(){gv_left = gv_right = 0; singular = false;}
        ~ThreeImplicitFunctions(){}
        
        virtual bool prepare_cell(int i, int j) = 0;

        virtual bool function_on_cell(double *val, int ir, int jr, int kl, int kr) = 0;

        GridValues * grid_value_left(void){return gv_left;}
        GridValues * grid_value_right(void){return gv_right;}

        bool is_singular(void){return singular;}
};

#endif // _THREEIMPLICITFUNCTIONS_


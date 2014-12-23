#ifndef _TWOIMPLICITFUNCTIONS_
#define _TWOIMPLICITFUNCTIONS_

#ifndef _INVALID_SEGMENT_
#define _INVALID_SEGMENT_ 0
#endif

#ifndef _VALID_SEGMENT_
#define _VALID_SEGMENT_ 1
#endif

#include <vector>
#include "GridValues.h"

// This is the template for any method which tries to use the Contour2p5.
// (At least) Three functions must be created for such end:
//
//     valid_segment(int i)     With i segment number (notice that Contour2p5 calls the point 2*i, thus
//                              here is used implicitly i as this first point and i+1 as the second point.
//     function_on_vertices(double *foncub, int i, int j, int kl)
//                              foncub is the array to be filled in the domain coordinates (i, j) given
//                              by the GridValues and kl is the index of the current segment, 0 and 1.
//     curve(...)               Is the method which be called for outsiders in order to evaluate the
//                              method. The signature have OPEN FORMAT and should be determine by the
//                              programer.

class TwoImplicitFunctions {
    private:
    protected:
        GridValues *gv;
        std::vector<RealVector> *oc; // Original curve
        bool singular;
    public:
        TwoImplicitFunctions(){gv = 0; singular = false;}
        ~TwoImplicitFunctions(){}
        
        virtual bool valid_segment(int i) = 0;

        virtual int function_on_vertices(double *foncub, int i, int j, int kl) = 0;


        GridValues * grid_value(void){return gv;}
        std::vector<RealVector> * segment_value(void){return oc;}
        
        bool is_singular(void){return singular;}
};

#endif // _TWOIMPLICITFUNCTIONS_


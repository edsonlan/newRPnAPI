#include "Quad3FluxParams.h"


//const double Quad3FluxParams::DEFAULT_A[3] = {0, 0 , 0};
//const double Quad3FluxParams::DEFAULT_B[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
//const double Quad3FluxParams::DEFAULT_C[3][3][3] = {
//                                                    {
//                                                     {0, 0, 0}, {0, 0, 0}, {0, 0, 0}
//                                                    },
//                                                    {
//                                                     {0, 0, 0}, {0, 0, 0}, {0, 0, 0}
//                                                    },
//                                                    {
//                                                     {0, 0, 0}, {0, 0, 0}, {0, 0, 0}
//                                                    }
//                                                   };
//

//Quad3FluxParams::Quad3FluxParams(void) : FluxParams(RealVector(27))
Quad3FluxParams::Quad3FluxParams(void) : FluxParams(RealVector(27))
{
    for (int i = 0; i < 27; i++) component(i, 0);
    component(0,3);
    component(2,1);
    component(10,1);

}

Quad3FluxParams::~Quad3FluxParams(){}

Quad3FluxParams::Quad3FluxParams(const RealVector & params) :FluxParams(params) 
{
    for (int i = 0; i < 27; i++) component(i, params(i));
}


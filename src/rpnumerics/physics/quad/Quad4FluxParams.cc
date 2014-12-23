#include "Quad4FluxParams.h"




const double Quad4FluxParams::DEFAULT_A[2] = { 0., 0. };
const double Quad4FluxParams::DEFAULT_B[2][2] = { { 0., .1 }, { -.1, 0. } };
const double Quad4FluxParams::DEFAULT_C[2][2][2] = { { { -1., 0. }, { 0., 1. } }, { { 0., 1. }, { 1., 0. } } };
//

Quad4FluxParams::Quad4FluxParams(void) :FluxParams(RealVector(56)) //TODO Criar os parametros (RealVector de muitas posicoes !!) de acordo com o que esta definido na fluxFunctionn
{
 
//    int m=2;

        for (int i = 0; i < 55; i++) {
        component(i, 1);
        }
//            component(i, DEFAULT_A[i]);
//            for (int j = 0; j < m; j++) {
//            component(m + i * m + j, DEFAULT_B[i] [j]);
//                for (int k = 0; k < m; k++) {
//                    component(m + m * m + i * m * m + j * m + k, DEFAULT_C[i] [j] [k]);
//                }
//            }
//        }
    
    
    
}

Quad4FluxParams::~Quad4FluxParams(){}

Quad4FluxParams::Quad4FluxParams(const RealVector & params) :FluxParams(params) 
{
}




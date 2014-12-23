#ifndef _VISCOSITY_MATRIX_
#define _VISCOSITY_MATRIX_

#include "ViscosityJetMatrix.h"
#include "RealVector.h"
#include "Matrix.h"

class Viscosity_Matrix {
    private:
    protected:
    public:
        Viscosity_Matrix(){}
        virtual ~Viscosity_Matrix(){}

        virtual void fill_viscous_matrix(const RealVector &p, ViscosityJetMatrix &m) const;
        virtual void fill_viscous_matrix(const RealVector &p, ViscosityJetMatrix &m, int degree) const;

        virtual bool is_constant(void){return true;}
        virtual bool is_identity(void){return true;}
};

#endif // _VISCOSITY_MATRIX_


#ifndef _VISCOSITYJETMATRIX_
#define _VISCOSITYJETMATRIX_

#include <vector>
#include "DoubleMatrix.h"

// Class to be used by Viscosity_Matrix as storage.
//
class ViscosityJetMatrix {
    private:
    protected:
        int nrows_, ncols_, nvar_;
    
        DoubleMatrix M_;

        //Matrix<DoubleMatrix > JM_;
        
        //std::vector<Matrix<DoubleMatrix > > HM_;

        void init();
    public:
        ViscosityJetMatrix();
        ViscosityJetMatrix(int nrows, int ncols, int nvar);
        ~ViscosityJetMatrix();

        DoubleMatrix       & M(void)       {return M_;}
        DoubleMatrix const & M(void) const {return M_;}

//        Matrix<DoubleMatrix >       & JM(void)       {return JM_;}
//        Matrix<DoubleMatrix > const & JM(void) const {return JM_;}

//        std::vector<Matrix<DoubleMatrix > >       & HM(void)       {return HM_;}
//        std::vector<Matrix<DoubleMatrix > > const & HM(void) const {return HM_;}
};

#endif // _VISCOSITYJETMATRIX_


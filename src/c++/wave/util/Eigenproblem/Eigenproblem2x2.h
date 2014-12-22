#ifndef _EIGENPROBLEM2X2_
#define _EIGENPROBLEM2X2_

#include "Eigenproblem.h"

class Eigenproblem2x2: public Eigenproblem {
    private:
    protected:
    public:
        Eigenproblem2x2();
        virtual ~Eigenproblem2x2();

        // Standard problem.
        //
        virtual void find_eigenpair(const DoubleMatrix &A, int index, Eigenpair &ep);
        virtual void find_eigenpairs(const DoubleMatrix &A, std::vector<Eigenpair> &eps);

        virtual void find_eigenvalue(const DoubleMatrix &A, int index, Eigenvalue &ev);
        virtual void find_eigenvalues(const DoubleMatrix &A, std::vector<Eigenvalue> &evs);

        // Generalized problem.
        //
        virtual void find_eigenpair(const DoubleMatrix &A, const DoubleMatrix &B, int index, Eigenpair &ep);
        virtual void find_eigenpairs(const DoubleMatrix &A, const DoubleMatrix &B, std::vector<Eigenpair> &eps);

        virtual void find_eigenvalue(const DoubleMatrix &A, const DoubleMatrix &B, int index, Eigenvalue &ev);
        virtual void find_eigenvalues(const DoubleMatrix &A, const DoubleMatrix &B, std::vector<Eigenvalue> &evs);
};

#endif // _EIGENPROBLEM2X2_


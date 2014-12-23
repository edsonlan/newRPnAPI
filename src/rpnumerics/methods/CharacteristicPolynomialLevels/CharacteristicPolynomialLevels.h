#ifndef _CHARACTERISTICPOLYNOMIALLEVELS_
#define _CHARACTERISTICPOLYNOMIALLEVELS_

#include "ImplicitFunction.h"
#include "ContourMethod.h"
#include "Coincidence.h"

#define LAMBDA_S 0
#define LAMBDA_E 1

// TODO: Separate the eigenvalue-related part.

class CharacteristicPolynomialLevels : public ImplicitFunction {
    private:
    protected:
        // Discriminat-related.
        //
        double discriminant_level;

        double discriminant(const DoubleMatrix &FJ, const DoubleMatrix &GJ);

        static int discriminant_function_on_square(CharacteristicPolynomialLevels *obj, double *foncub, int i, int j);

        // Discriminat derivative-related.
        //
        int index_derivative_discriminant_variable;

        RealVector abc(const DoubleMatrix &FJ, const DoubleMatrix &GJ);
        double localm_(const RealVector &xa, const RealVector &xb, double flip, double t, RealVector &xn);

        static double max_discriminant(CharacteristicPolynomialLevels *obj, int i, int j);
        static int derivative_discriminant_function_on_square(CharacteristicPolynomialLevels *obj, double *foncub, int i, int j);

        // Eigenvalue-related.
        //
        double eigenvalue_level;
        int eigenvalue_family;
        DoubleMatrix eigenvalue_matrix;

        static int eigenvalue_function_on_square(CharacteristicPolynomialLevels *obj, double *foncub, int i, int j);

        // Common.
        //
        const FluxFunction *ff;
        const AccumulationFunction *aa;
        bool is_improvable;

        int (*fos)(CharacteristicPolynomialLevels *obj, double *foncub, int i, int j);
    public:
        // Discriminat-related.
        //
        double discriminant(const FluxFunction *f, const AccumulationFunction *a, const RealVector &p);

        int discriminant_curve(const FluxFunction *f, const AccumulationFunction *a, 
                               GridValues &g, 
                               double level,
                               std::vector<RealVector> &curve);

        void discriminant_curve(const FluxFunction *f, const AccumulationFunction *a, 
                                GridValues &g, 
                                const std::vector<double> &level,
                                std::vector<std::vector<RealVector> > &curve);

        int discriminant_curve(const FluxFunction *f, const AccumulationFunction *a, 
                               GridValues &g, 
                               const RealVector &p,
                               std::vector<RealVector> &curve, double &lev);

        // Discriminat derivative-related.
        //

        // TODO:
        // The value of u should be selected internally and not by the user.
        //
        int derivative_discriminant_curve(const FluxFunction *f, const AccumulationFunction *a, 
                                          GridValues &g, int u,
                                          std::vector<RealVector> &curve, 
                                          std::vector<double> &discriminant_on_segment);

        // Derivative of the discriminant with respect to the variable indexed by u.
        //
        double derivative_discriminant(const FluxFunction *f, const AccumulationFunction *a, const RealVector &p, int u);

        int complete(const RealVector &p0, const RealVector &p1, const RealVector &p_init, RealVector &p_completed);

        // Eigenvalue-related.
        //
        double eigenvalue(const FluxFunction *f, const AccumulationFunction *a,
                          const RealVector &p, int fam);

        int eigenvalue_curve(const FluxFunction *f, const AccumulationFunction *a, 
                             GridValues &g, 
                             double level, int family, 
                             std::vector<RealVector> &curve);

        // With Coincidence
        int eigenvalue_curve(const Coincidence *c,
                             const FluxFunction *f, const AccumulationFunction *a, 
                             GridValues &g, 
                             const RealVector &p, int family, 
                             std::vector<RealVector> &curve);
        static int exact_eigenvalue_function_on_square(CharacteristicPolynomialLevels *obj, double *foncub, int i, int j);

        void eigenvalue_curve(const FluxFunction *f, const AccumulationFunction *a, 
                              GridValues &g, 
                              const std::vector<double> &level, int family, 
                              std::vector<std::vector<RealVector> > &curve);

        int eigenvalue_curve(const FluxFunction *f, const AccumulationFunction *a, 
                             GridValues &g, 
                             const RealVector &p, int family, 
                             std::vector<RealVector> &curve, double &lev);

        // Common.
        //
        CharacteristicPolynomialLevels();
        virtual ~CharacteristicPolynomialLevels();

        int function_on_square(double *foncub, int i, int j){return (*fos)(this, foncub, i, j);}
        bool improvable(void){return is_improvable;}
};

#endif // _CHARACTERISTICPOLYNOMIALLEVELS_


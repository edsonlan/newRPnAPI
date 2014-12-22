#ifndef _EIGEN_
#define _EIGEN_

#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "RealVector.h"
#include "DoubleMatrix.h"

// Forward declarations
class FluxFunction;
class AccumulationFunction;

using namespace std;

// When beta is near zero (for the generalized eigenproblem),
// function eig() must return this value. This situation is anomalous.
//
#define _EIGEN_BETA_NEAR_ZERO_ -123
#define COMPLEX_EIGENVALUE (-7)
#define SUCCESSFUL_PROCEDURE 2
#define ABORTED_PROCEDURE (-7)
#define LAMBDA_ERROR (-7)
#define LAMBDA_NOT_INCREASING (-7)
#define LAMBDA_NOT_DECREASING (-7)
#define _SIMPLE_ACCUMULATION_  10  // Traditional rarefaction, using dgeev.
#define _GENERAL_ACCUMULATION_ 11  // Rarefactio


// Eigenproblem, A*v = lambda*v
//
extern "C" void dgeev_(const char*, const char*, int*, double*, int*, double*, double*, 
           double*, int*, double*, int*, double*, int*, 
           int*);

// Generalized eigenproblem, A*v = lambda*B*v
//
extern "C" void dggev_(const char*, const char*,  // JOBVL, JOBVR
                  int*,                      // N
                  double*, int*,             // A, LDA
                  double*, int*,             // B, LDB
                  double*,                   // ALPHAR
                  double*,                   // ALPHAI
                  double*,                   // BETA
                  double*, int*,             // VL, LDVL,
                  double*, int*,             // VR, LDVR,
                  double*, int*,             // WORK, LWORK
                  int*                       // INFO
                 );

/* Struct to hold an eigenpair. */
struct eigenpair {
    public:
        eigenpair(){
        }

        eigenpair(int n){
            r = 0.0;
            i = 0.0;

            vlr.resize(n);
            vli.resize(n);
            vrr.resize(n);
            vri.resize(n);

            for (int i = 0; i < n; i++){
                vlr[i] = vli[i] = vrr[i] = vri[i] = 0.0;
            }
        }

        double r;           // Real part of the eigenvalue
        double i;           // Imaginary part of the eigenvalue

        vector<double> vlr; // Real part of the left-eigenvector
        vector<double> vli; // Imaginary part of the left-eigenvector
        vector<double> vrr; // Real part of the right-eigenvector
        vector<double> vri; // Imaginary part of the right-eigenvector
        
       eigenpair operator=(const eigenpair &original){
           if (this != &original){
               r = original.r;
               i = original.i;
            
               int n = original.vlr.size();
               vlr.resize(n);
               vli.resize(n);
               vrr.resize(n);
               vri.resize(n);
            
               for (int ii = 0; ii < n; ii++){
                   vlr[ii] = original.vlr[ii];
                   vli[ii] = original.vli[ii];
                   vrr[ii] = original.vrr[ii];
                   vri[ii] = original.vri[ii];
               }
            }

            return *this;
        }
        
};
/* Struct to hold an eigenpair. */

/* Class Eigen. */
class Eigen {
    private:
        static double epsilon;
        static void transpose(int, double*);
        static void fill_eigen(int, eigenpair*, double*, double*, double*, double*);
        static int eigen_comp(const void *, const void *); // To be deprecated
        static bool eigen_compare(const eigenpair&, const eigenpair&);
        static void sort_eigen(int, eigenpair*);


        static bool eigen_compare_using_eigenvalues(const eigenpair&, const eigenpair&);
        static bool eigen_compare_using_eigenvectors(const eigenpair&, const eigenpair&);

        static bool (*eigen_sort_function)(const eigenpair&, const eigenpair&);
    protected:
    public:
        static int eig(int n, const double*, vector<eigenpair>&);                // Eigenproblem
        static int eig(int n, const double*, const double*, vector<eigenpair>&); // Generalized eigenproblem
        static int eig(int n, const double*, const double*, int family, eigenpair&); // Generalized eigenproblem for a given family
        static int eig(int n, const double *A, const double *B, int family, RealVector &r); //Generalized eigenproblem for a given family, returning ONLY the right-eigenvector
        static int eig(const RealVector &p, const FluxFunction *f, const AccumulationFunction *g, std::vector<eigenpair> &e);

        static void print_eigen(const vector<eigenpair> &);
        static void eps(double);
        static double eps(void);

        static void fill_eigenpairs(const FluxFunction *f, const AccumulationFunction *a, const RealVector &u, std::vector<eigenpair> &e);
        static void fill_eigenvalues(const FluxFunction *f, const AccumulationFunction *a, const RealVector &u, std::vector<double> &lambda);

        // Get/set method to order a list of eigenpairs.
        //
        static void list_order_eigenpairs(std::vector<bool (*)(const eigenpair&, const eigenpair&)> &order_function, std::vector<std::string> &name){
            order_function.clear();
            name.clear();

            order_function.push_back(&eigen_compare_using_eigenvalues);
            name.push_back(std::string("Using families"));

            order_function.push_back(&eigen_compare_using_eigenvectors);
            name.push_back(std::string("Using types"));
        }

        static void set_order_eigenpairs(bool (*f)(const eigenpair&, const eigenpair&)){
            eigen_sort_function = f;

            return;
        }

};
/* Class Eigen. */


#endif // _EIGEN_

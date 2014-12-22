#ifndef _GRIDVALUES_
#define _GRIDVALUES_

#include <stdio.h>
#include "FluxFunction.h"
#include "AccumulationFunction.h"
#include "RealVector.h"
#include "eigen.h"
//#include "Boundary.h"

class Boundary;

#include "Matrix.h"
#include "DoubleMatrix.h"

// The definitions below will help the Contour to decide whether a cell is to be
// processed (and if so, how) or not.
//
#ifndef CELL_IS_INVALID
#define CELL_IS_INVALID 0
#endif

#ifndef CELL_IS_TRIANGLE    // Lower triangle
#define CELL_IS_TRIANGLE 1
#endif

#ifndef CELL_IS_SQUARE
#define CELL_IS_SQUARE 2
#endif

#define UP    10
#define LEFT  11
#define DOWN  12
#define RIGHT 13
#define NONE  14

class GridValues {
    private:
    protected:
        // Fill the bare minimum of the grid values.
        //
        void fill_values_on_grid(const Boundary *b, 
                                 const RealVector &min, const RealVector &max,
                                 const std::vector<int> &number_of_cells);

        // For a double-linked list.
        //
        static GridValues *first_, *last_, *active_;

        GridValues *prev_, *next_;
    public:
        // Values on the grid. Every value MUST have a boolean companion that can be used
        // to know if the values are already computed or not.

        Matrix<RealVector>               grid;                       // Grid proper
        Matrix<bool>                     point_inside;               // The point belongs to the domain (verified via Boundary)
        Matrix<int>                      cell_type;                  // One of CELL_IS_INVALID, CELL_IS_TRIANGLE (x)or CELL_IS_SQUARE.
        bool                             grid_computed;              // Already computed?

        Matrix<RealVector>               F_on_grid;                  // Flux
        Matrix<RealVector>               G_on_grid;                  // Accumulation
        bool                             functions_on_grid_computed; // Already computed?

//        Matrix< Matrix<double> >         JF_on_grid;                 // Jacobians of the flux
//        Matrix< Matrix<double> >         JG_on_grid;                 // Jacobians of the accumulation

        Matrix<DoubleMatrix>             JF_on_grid;                 // Jacobians of the flux
        Matrix<DoubleMatrix>             JG_on_grid;                 // Jacobians of the accumulation
        bool                             Jacobians_on_grid_computed; // Already computed?

        Matrix< std::vector<double> >    dd;                         // Directional derivatives
        bool                             dd_computed;                // Already computed?

        Matrix< std::vector<eigenpair> > e;                          // Eigenpairs
        Matrix< std::vector<bool> >      eig_is_real;                // Are the eigenvalue real or complex?
        Matrix<bool>                     cell_is_real;               // Is the whole cell real or complex?
        bool                             e_computed;                 // Already computed?

        void clear_computations(){
            functions_on_grid_computed = false;
            Jacobians_on_grid_computed = false;
            dd_computed = false;
            e_computed = false;

            return;
        }

        // TODO:
        // Replace (or complement) the use of cell_is_real by (with) a Matrix of discriminants.

        RealVector grid_resolution;                                  //Number of cells

        // Here the future programmers will add the values that will prove necessary later.

        // Constructor. Fills the spatial grid.
        //
        GridValues(const Boundary *b, 
                   const RealVector &pmin, const RealVector &pmax,
                   const std::vector<int> &number_of_cells);

        ~GridValues();

        // Set the grid. This function could be merged with fill_values_on_grid().
        //
        void set_grid(const Boundary *b, 
                      const RealVector &min, const RealVector &max,
                      const std::vector<int> &number_of_cells);

        // Fill F and G.
        // 
        void fill_functions_on_grid(const FluxFunction *ff, const AccumulationFunction *aa);

        // Fill JF and JG.
        // 
        void fill_Jacobians_on_grid(const FluxFunction *ff, const AccumulationFunction *aa);

        // Fill eigenpairs.
        //
        void fill_eigenpairs_on_grid(const FluxFunction *ff, const AccumulationFunction *aa);

        // Fill the directional derivatives.
        //
        void fill_dirdrv_on_grid(const FluxFunction *ff, const AccumulationFunction *aa);

        bool inside(const RealVector &p);

        bool cell(const RealVector &p, int &row, int &col);

        // Set/get the active object.
        //
        static void active(GridValues *g);
        static GridValues* active();

        // First/last object.
        //
        static GridValues* first();
        static GridValues* last();

        // Access to the object's next and previous objects.
        //
        GridValues* next();
        GridValues* prev();

        const GridValues* next() const;
        const GridValues* prev() const;
};

#endif // _GRIDVALUES_


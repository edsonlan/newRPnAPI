/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) HyperCube.h
 */

#ifndef _HyperCube_H
#define _HyperCube_H
#define ERROR -666
#define OK 0
#define KLUDGE 1e-6
/*
 * ---------------------------------------------------------------
 * Includes:
 */

#include <stdio.h>
#include <math.h>
#include "eigen.h"
#include <iostream>

#include "Matrix.h"

/*
 * ---------------------------------------------------------------
 * Definitions:
 */
extern "C" {


    void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);


}

class HyperCube {
private:




    /*
     * Instituto de Matematica Pura e Aplicada - IMPA
     * Departamento de Dinamica dos Fluidos
     *
     */


    /*
     * If n is out of range, factorial returns 0.
     */
    int factorial(int n);


    int pow2(int n);


    void Eval_dimf(int *dimFACE);




    int CubeFace(int dim, int faceDim /*Provalmente faltam ponteros*/);

    /*
     * @param altFace The new face array
     * @throws IllegalArgumentException If the array has an invalid size.
     */
    /*  TODO Parece ser que nao temos como fazer essa mensagem de erro funcionar.
    void setFaceArray(int[] [] altFace) {
        //Check sizes

    /*
     * This routine creates the standard hypercube vertices (in cvert_),
     * the basic simplex vertices (in bsvert_), and the permutations of 1, 2, ..., n_  (in permutations).
     */



    /*
     * This routine creates the array  perm_  containing the permutations
     * of the numbers from 1 to  n_.  It also sets  nsimp_  to be  n_  factorial.
     */
    void mkperm(int *perm_, // perm_[n_][nsimp_]
            int n_, int nsimp_);



    /** This routine makes the arrays  face_  and  facptr_. It also defines the variable  nface_. */
    int mkfcfp(int *face_, int *facptr_, // face_[m_ + 1][dimf_], facptr_[nsimp_][nsface_]
            int nsface_, int dimf_, int nsimp_, int m_, int n_,
            int *bsvert_, int *comb_, // bsvert_[n_ + 1][n_], comb_[numberOfCombinations][m_ + 1]
            int *perm_, int *stor); // perm_[n_][nsimp_], storm[m_ + 1]





    /**
     * This routine searches the columns of a2  for the column array a1.
     * @param a1 The column to be searched
     * @param a2 The array to be analized
     * @param n1 Number of columns of the array
     * @param n2 Number of lines of the array
     * @returns  <p> -1  ----------> a1  not found in  a2 </p> <p> column no.  --> a1  found in  a2 </p>
     */

    int search(int *a1, int *a2, int n1, int n2, int n3);


    /*
     * This routine makes the adjacency array  fnbr_  for the faces in
     * the basic simplex. Two faces in comb_ are adjacent if they differ in exactly one vertex.
     */
    void mkfnbr(int *fnbr_, int *comb_, //fnbr_[nsface_][nsface_], comb_[numberOfCombinations][m_ + 1]
            int n_, int m_, int nsface_,
            int *stor); //storn[n_+1]


    /**
     * This function determines which face_ is opposite face_  f  in
     * a cube.  It is assumed that face_  f  lies on a cube face_.
     * @param f Face number
     * @param xind Indicates which component is constant.
     * @param value Is that constant value.
     * @param oppvrt Auxiliar vector
     * @returns See the search method return values
     */
    int mkoppf(int f, int xind, int value, int n_, int m_, int nface_, int dimf_, int *face_, int *oppvrt);


    /**
     * Funcao para Debug!! Imprime array de inteiros trocando linhas com colunas.
     * <p><h3> Deve ser eliminada apos testes da classe </h3></p>
     * @param name Variable name
     * @param matrix The bidimensional array to be printed
     * @param k1 Size of the array
     * @param k2 Size of the array
     */
    /*void printInt(String name, int[] [] matrix, int k1, int k2) {
        int i, j;
        System.out.println(name + ": ");
        for (i = 0; i < k1; i++) {
            for (j = 0; j < k2; j++) {
                System.out.print(matrix[i] [j] + " ");
            }
            System.out.println("");
        }
    }
     */
    /**
     * Funcao para Debug!! Imprime array de doubles na ordem correta.
     * <p><h3> Deve ser eliminada apos testes da classe </h3></p>
     * @param name Variable name
     * @param matrix The bidimensional array to be printed
     * @param k1 Size of the array
     * @param k2 Size of the array
     */
    /*void printDouble(String name, double[] [] matrix, int k1, int k2) {
        int i, j;
        System.out.println(name + ": ");
        for (i = 0; i < k1; i++) {
            for (j = 0; j < k2; j++) {
                System.out.print(matrix[i] [j] + " ");
            }
            System.out.println("");
        }
    }
     */
    /**
     * Doubles the size of an array
     * @param source Array de origem
     * @return Returns an array containing the elements of source with twice its size
     */
    /*Object doubleArray(Object source) {
        int sourceLength = Array.getLength(source);
        Class arrayClass = source.getClass();
        Class componentClass = arrayClass.getComponentType();
        Object result = Array.newInstance(componentClass, sourceLength * 2);
        System.arraycopy(source, 0, result, 0, sourceLength);
        return result;
    }
     */




    /**
     * This routine calls  affslv  to obtain the solutions on the faces in  face.
     * @param sol_ Array of solutions
     * @param exstfc Bit array indicating which faces in array "face" are really to be considered by cubsol.
     * Dimension of exstfc should be nface
     */
    int mksoln(double *sol_, int dims_, int *sptr_, int nsoln_, double *foncub, int *exstfc, int *face, int dimf_, double *cvert, int ncvert_, int n_, int m_, int nface, double *u, double *g, double *x, int *wrki);


    //    double auxsol[n_]; NAO PARECE SER USADA!!

    /*	    int wrki[m_];
                double u[n_][m_ + 1];
                double g[m_][m_];
                double gx[m_];
                double x[m_];
                double solution[n_][m_];
     */

    //ATENCAO em Fortran eram I/O
    //    int nface = cFace_.getNumberOfFaces();
    //    int[] [] face;
    //    face = cFace_.getFaceArray();
    //    double[] [] cvert;
    //TODO em Fortran era I/O    cvert = cFace_.getVertexCoordinatesArray();
    //FECHA ATENCAO

    //loop over the faces to find each solution (if there is one)
    /*java:    nsoln_ = 0;
        for (indf = 0; indf < nface; indf++) {
            sptr_[indf] = 0;
            //if the indf-face is not to be considered by mksoln, skip it
            if (exstfc[indf] != 0) {
                    //boolean zeroFound = false;
                //initialize gx
                v = face[0*dimf_ + indf];
                for (i = 0; i < m_; i++) {
                    gx[i] = foncub[i*m_ + v];
                }

                //set the function values  at the face vertices
                for (k = 0; k < m_; k++) {
                    v = face[(k + 1)*dimf_ + indf];
                    for (i = 0; i < m_; i++) {
                        g[i* m_ +k] = foncub[i*m_ + v + 1]; //TODO: Porque nao falta um +1 apos v?
                    }
                }

                //if (zeroFound) {
                //	zeroCount++;
                //}
                //call the affine solver to solve in the standard simplex
                //if(!zeroFound || (zeroCount <= 1)) {
                            flag = affslv(x, gx, g, m_, wrki);
                            //skip the rest if no solution ( note -- pointer initialized to 0 )

                            // solve the system
                            // flag = calculateSolution(x, gx, g, m_, wrki);

                            if (flag == 0) {
                                //set the pointer to the solution
                                nsoln_ = nsoln_ + 1;
                                if (dims_ < nsoln_) {
                                    // TODO: Originally in Fortran mksoln returns 1 as an error (?).
                                    return 1;
                                    //TODO.... Uma tabela vai ser precisa para isto.... (?)
                                    //         Nao ha como saber qual o valor inicial de dims_
    //	                        dims_ *= 2;
    //	                        for (j = 0; j < n_; j++) {
    //	                            sol_[j] = (double[]) doubleArray(sol_[j]);

    //	                        }
    //

    // TODO: Parece faltar sptr_[indf] = nsoln_;

                                }
                                //set the face vertices
                                for (k = 0; k < m_ + 1; k++) {
                                    v = face[k*dimf_ + indf];
                                    for (i = 0; i < n_; i++) {
                                        u[i*(m_+1) + k] = cvert[(v + 1)*n_ + i];//TODO: Prestar atencao com a ordem a direita, no cvert!!
                                    }
                                }
                                //transform the solution back from the standard simplex
                                double tmp[n_];
                                tmp[0]=-10;
                                afftrn(&tmp[0], u, x, n_, m_);

                                if(tmp[0] != -10) {
                                    sptr_[indf] = nsoln_;
                                    int pont_solution;
                                        for(pont_solution = 0; pont_solution < n_; pont_solution++) {
                                            sol_[pont_solution*dims_ + (nsoln_ - 1)] = tmp[pont_solution];
                                        }
                                } else {
                                    // TODO: In the same vein, we return 1 as error.
                                    return 1;
                                }
                            }
                //}

            }
        }
        return 0;
    }java:*/

    /**
     * This routine finds the zero  x  of the affine function whose values at
     * the vertices of the standard m_-simplex are stored in the columns of  g.
     * It also checks whether the solution is inside the simplex and sets flag appropriately.
     * @param x
     * @param g0
     * @param g
     * @param wrki
     * @return  <ul> <p> -1      Error in solving the system </p> <p>  0      Valid solution inside simplex </p>
     * <p>  1      Valid solution outside simplex </p> <ul>
     */
    /* void affslv(double[] x, double[] g0, double[] [] g, int[] wrki) {
     //local variables
     int i, j;
     //subtract  g0  from  g1, g2, ..., gm
     for (j = 0; j < m_; j++) {
         for (i = 0; i < m_; i++) {
             g[i] [j] = g[i] [j] - g0[i];
         }
     }
     //change sign of  g0
     for (i = 0; i < m_; i++) {
         g0[i] = -g0[i];
     }

 }*/

    int affslv(double *x, double *g, int m_, int *wrki);



    /**
     * <p> This routine maps  x  in the standard m_-simplex to  y  in the
     * n_-dimensional m_-simplex with vertices  u(.,1), u(.,2), ...,  u(.,m_+1). </p>
     * <p> The map is the affine transformation which takes the vertices of the
     * standard simplex to those contained in the columns of  u. </p>
     * @param y
     * @param u
     * @param x
     */
    void afftrn(double *solution, double *u, double *x, int n_, int m_);



    int solver(int n, double *A, double *b, double *x);





    void smpptr(int *solptr_, int *sptr_, int *facptr, int nsimp_, int nsface_);




    /**
         the following routine builds the array of pointers to the edges in
    c     the solution polyhedron which are parallel to a coordinate plane.
    c
    c   input
    c   -----
    c     solptr  --  a 2-dimensional array containing pointers to the
    c                 solution in  sol.  i.e.,  solptr(i,j)  is the
    c                 column of  sol  containing the solution for the
    c                 i-th face in the j-th simplex.  (if  solptr  is
    c                 0, the corresponding face contains no solution.)
    c
    c     fnbr    --  a 2-dimensional array containing the adjacency
    c                 information for the faces in the basic simplex.
    c                 i.e.,  fnbr(i,j)  is 1 (resp. 0) if the i-th and
    c                 j-th faces are (resp. are not) adjacent in the
    c                 basic simplex.  here,  i,j  refer to the order
    c                 determined by the array  comb.
    c
    c     nsimp   --  the number of simplices in the hypercube.
    c
    c     nsface  --  the number of faces in the basic simplex.
    c
    c   output
    c   ------
    c     edges   --  a 2-dimensional array each column of which contains
    c                 pointers to the array  sol  at which the vertices
    c                 of an edge of the solution polyhedron are stored.
    c
    c     nedges  --  the number of edges found.
    c
    c     smpedg  --  a 2-dimensional array each column of which contains
    c                 the beginning and ending indices for the array
    c                 edges.  i.e., the edges for simplex number  ns
    c                 are in columns  smpedg(1,ns)  to  smpedg(2,ns)
    c                 of  edges.  there are no edges if
    c                        smpedg(1,ns) .gt. smpedg(2,ns).
    c
     */
    int mklevl(int *edges_, int dime_, int nedges_, int *smpedg_, int *solptr_, int *fnbr_, int nsimp, int nsface, int *facatt, int *ptfatt, int *facptr, int fdim_, int nface_);


    int mkpoly(int *polyv_, int dimv, int npolys, int *smppol_, int *solptr_, int *fnbr_, int nsimp_, int nsface_, int dimp_);


    // The following routine prints an integer matrix, priorizing the longest lenght as lines.
    // For traking the actual size is printed within the name

    //    void putmi(const char *name, int *matrix_, int k1_, int k2_);
    //
    //
    //    void putmi2(const char *name, int *matrix_, int k1_, int k2_, int kmax);

    // The following routine prints a real matrix, priorizing the longest lenght as lines.
    // For traking the actual size is printed within the name

    //    void putmf(const char *name, double *matrix_, int k1_, int k2_);


    //    void putmf2(const char *name, double *matrix_, int k1_, int k2_, int kmax);


    //    void putme(const char *name, double *matrix_, int k1_, int k2_);




public:
    void cpp_mkcube(Matrix<double> &cpp_cvert_, // cvert_[ncvert_][n_]
            Matrix<int> &cpp_bsvert_, // bsvert_[n_ + 1][n_]
            Matrix<int> &cpp_perm_, // perm_[n_][nsimp_]
            int ncvert_, int nsimp_, int n_);

    void mkcube(double *cvert_, // cvert_[ncvert_][n_]
            int *bsvert_, // bsvert_[n_ + 1][n_]
            int *perm_, // perm_[n_][nsimp_]
            int ncvert_, int nsimp_, int n_);






    /*
     * This routine creates the list  face_  of all m_-dimensional faces
     * of the  n_!  simplices obtained by permuting the coordinates of
     * the basic simplex.  to do this it needs the array  comb_  of all
     * subsets of size  m_+1  from a set of size  n_+1.  since  face_  contains
     * only distinct faces, a pointer array  facptr_  is needed to
     * determine which column of  face_  corresponds to each permutation and combination.
     */

    int mkface(int *face_, int *facptr_, int *fnbr_, // face_[m_ + 1][dimf_], facptr_[nsimp_][nsface_]
            int dimf_, int nsimp_, int n_, int m_, int nsface_,
            int *bsvert_, int *comb_, // bsvert_[n_ + 1][n_], comb_[numberOfCombinations][m_ + 1]
            int *perm_, int *storn_, int *storm_); // perm[], storn[n_ + 1], storm[m_ + 1]






    /*AQUI COMECA A PARTE CORRESPONDENTE A CUBSOLVER*/

    /*
     * Instituto de Matematica Pura e Aplicada - IMPA
     * Departamento de Dinamica dos Fluidos
     *
     */

    // TODO: exstfc eh usado dentro de mksoln, porem ainda nao foi inicializado
    int cpp_cubsol(int *solptr_, Matrix<double> &cpp_sol_, int dims_, int *sptr_, int nsoln_,
            double *foncub, int *exstfc, int *face,
            int *facptr_, int dimf_, double *cvert, int ncvert_, int n_, int m_,
            int nsimp_, int nsface_, int nface, double *u, double *g,
            double *x, int *wrki);

    int cubsol(int *solptr_, double *sol_, int dims_, int *sptr_, int nsoln_,
            double *foncub, int *exstfc, int *face,
            int *facptr_, int dimf_, double *cvert, int ncvert_, int n_, int m_,
            int nsimp_, int nsface_, int nface, double *u, double *g,
            double *x, int *wrki);





    /** The following routine builds the array of pointers to the edges_ in the solution polyhedron. */

    /*input
       -----
         solptr  --  a 2-dimensional array containing pointers to the
                     solution in  sol.  i.e.,  solptr(i,j)  is the
                     column of  sol  containing the solution for the
                     i-th face in the j-th simplex.  (if  solptr  is
                     0, the corresponding face contains no solution.)

         fnbr    --  a 2-dimensional array containing the adjacency
                     information for the faces in the basic simplex.
                     i.e.,  fnbr(i,j)  is 1 (resp. 0) if the i-th and
                     j-th faces are (resp. are not) adjacent in the
                     basic simplex.  here,  i,j  refer to the order
                     determined by the array  comb.

         nsimp   --  the number of simplices in the hypercube.

         nsface  --  the number of faces in the basic simplex.

       output
       ------
         edges   --  a 2-dimensional array each column of which contains
                     pointers to the array  sol  at which the vertices
                     of an edge of the solution polyhedron are stored.

         nedges  --  the number of edges found.

         smpedg  --  a 2-dimensional array each column of which contains
                     the beginning and ending indices for the array
                     edges.  i.e., the edges for simplex number  ns
                     are in columns  smpedg(1,ns)  to  smpedg(2,ns)
                     of  edges.  there are no edges if
                            smpedg(1,ns) .gt. smpedg(2,ns).
    c*/
    int combination(int n, int m);
    void findInds(int *faceptInd_, int *fdim_, //pointers to initializate these values
            int n_, int m_,
            double *cvert_, int *face_, // cvert_[ncvert_][n_], face_[m_ + 1][dimf_]
            int dimf_, int nface_, int nsimp_, int nsface_, int *facptr_);

    int mkcomb(int *comb_, int np, int mp);

    int cpp_mkedge(Matrix<int> &cpp_edges_, int dime_, int nedges_, int *smpedg_, int *solptr_, int *fnbr_, int nsimp_, int nsface_);

    int mkedge(int *edges_, int dime_, int nedges_, int *smpedg_, int *solptr_, int *fnbr_, int nsimp_, int nsface_);

    void mkflst(int *facept_, int *ptrf_, // facept_[faceptInd_][2], ptrf_[nface_][2]
            int *facatt_, int *ptfatt_, // facatt_[fdim_][2], ptfatt_[nface_][2]
            int dimf_, int fdim_, int n_, int m_, int nface_, int nsface_, int nsimp_, int ncvert_, int faceptInd_,
            double *cvert_, int *face_, // cvert_[ncvert_][n_], face_[m_ + 1][dimf_]
            int *facptr_, int *wrk_);



    void putmf(const char *name, double *matrix_, int k1_, int k2_);
    void putmf2(const char *name, double *matrix_, int k1_, int k2_, int kmax);
    void putmi(const char *name, int *matrix_, int k1_, int k2_);
    void putmi2(const char *name, int *matrix_, int k1_, int k2_, int kmax);

    HyperCube(){};
    ~HyperCube(){};
};

#endif //! _HyperCube_H

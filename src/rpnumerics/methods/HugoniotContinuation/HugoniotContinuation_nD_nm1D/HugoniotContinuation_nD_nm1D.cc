#include "HugoniotContinuation_nD_nm1D.h"

// Dan and Morante finished this code on Sept. 5, 2013. It should work ok when JetMatrix deals correctly
// with different number of equations and of variables.
//
//
// Equation: \frac{\partial}{\partial t} G(V) + \frac{\partial}{\partial x} u F(V) = 0,
// 
// where V \in R^{}
//
void HugoniotContinuation_nD_nm1D::jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                                                const RealVector &G, const DoubleMatrix &JG,
                                                const RealVector &C, const DoubleMatrix &JC,
                                                RealVector &H, DoubleMatrix &nablaH){
                                              
//void HugoniotContinuation_nD_nm1D::jet_Hugoniot(const JetMatrix &JM_F, const JetMatrix &JM_G, 
//                                                /*RealVector &H, DoubleMatrix &nablaH*/
//                                                DoubleMatrix &Hugoniot_matrix, 
//                                                std::vector<DoubleMatrix> &nabla_Hugoniot_matrix,
//                                                Matrix<DoubleMatrix> &nabla_nabla_Hugoniot_matrix){

    // The number of equations is the number of rows of the Jacobians.    
    int n = JF.rows();

    // Notice that in this case the dimension of the domain is one less than the dimension of the image.
    // Therefore, Jacobians should be of size: n, n - 1. JetMatrix deals with this kind
    // of situations.

    // [F] & [G]
    //
    RealVector diff_G = G - ref.G;
    RealVector diff_F = F - ref.F;
    
    // Eq. 2.1
    //
    DoubleMatrix Rankine_Hugoniot_equations(n, 3);
    for (int i = 0; i < n; i++){
        Rankine_Hugoniot_equations(i, 0) = diff_G(i);
        Rankine_Hugoniot_equations(i, 1) = -F(i); // <=======
        Rankine_Hugoniot_equations(i, 2) = ref.F(i);
    }
    
    // Nabla of Eq. 2.1
    //
    std::vector<DoubleMatrix> nabla_Rankine_Hugoniot_equations(n - 1);
    for (int k = 0; k < n - 1; k++){
        nabla_Rankine_Hugoniot_equations[k].resize(n, 3); // Possibly: resize(n, 2), since the last column is of zeroes.
        
        for (int i = 0; i < n; i++){
            nabla_Rankine_Hugoniot_equations[k](i, 0) = JG(i, k);
            nabla_Rankine_Hugoniot_equations[k](i, 1) = -JF(i, k);
            nabla_Rankine_Hugoniot_equations[k](i, 2) = 0.0;
        }
    }

    // Eqs. 2.3 and 2.4.
    //
    // Compute the determinants of the minors.
    //
    std::vector<DoubleMatrix> Rankine_Hugoniot_minor(n - 2);
    H.resize(n - 2);
    
    for (int i = 0; i < n - 2; i++){
        Rankine_Hugoniot_minor[i].resize(3, 3); // Possibly: resize(n, 2), since the last column is of zeroes.
        
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){  // Possibly: k < 2, since the last column is of zeroes.
                Rankine_Hugoniot_minor[i](j, k) = Rankine_Hugoniot_equations(i + j, k); // Eq. 2.3 proper.
            }
        }
        
        H(i) = det(Rankine_Hugoniot_minor[i]); // Eq. 2.4 proper.
    }

    // Eq. 2.5.
    //
    // The first index is related to the number of partial derivatives.
    // The second index is related to the number of minors.
    //
    Matrix<DoubleMatrix> nabla_minor_Hugoniot_matrix(n - 1, n - 2); 
    for (int i = 0; i < n - 1; i++){
        for (int m = 0; m < n - 2; m++){ // Minors
            nabla_minor_Hugoniot_matrix(i, m).resize(3, 3);  // Possibly: resize(n, 2), since the last column is of zeroes.
        
            for (int j = 0; j < 3; j++){
                for (int k = 0; k < 3; k++){
                    nabla_minor_Hugoniot_matrix(i, m)(j, k) = nabla_Rankine_Hugoniot_equations[i](i + j, k);
                }
            }
        }
    }

    // Eq. 2.6.
    //
    // nablaH.
    //
    nablaH.resize(n - 1, n - 2);
    for (int i = 0; i < n - 1; i++){
        for (int m = 0; m < n - 2; m++){
            //nabla_determinant_minor_Hugoniot_matrix(i, m) = derivative_det(minor_Hugoniot_matrix[m], nabla_minor_Hugoniot_matrix(i, m));
            nablaH(i, m) = derivative_det(Rankine_Hugoniot_minor[m], nabla_minor_Hugoniot_matrix(i, m));  // The new version of derivative_det would be used.
        }
    }

    //nabla_nabla_Hugoniot_matrix.resize(n - 1, n - 1);
    /*
    for (int m = 0; m < n - 1; m++){
        for (int p = 0; p < n - 1; p++){
            nabla_nabla_Hugoniot_matrix(m, p).resize(n, 3);
            
            DoubleMatrix J2G = JM_G.extract_matrix_from_Hessian(p);
            DoubleMatrix J2F = JM_F.extract_matrix_from_Hessian(p);
            
            for (int i = 0; i < n; i++){
                nabla_nabla_Hugoniot_matrix(m, p)(i, 0) = ;
                nabla_nabla_Hugoniot_matrix(m, p)(i, 1) = ;
                nabla_nabla_Hugoniot_matrix(m, p)(i, 2) = ;
            }
        }
    }
    */
/*


    Matrix<DoubleMatrix> derivative_minor_Hugoniot_matrix(n - 2, n - 1);
    for (int row = 0; row < n - 2; row++){
        for (int col = 0; col < n - 1; col++){
            derivative_minor_Hugoniot_matrix(row, col).resize(3, 3);
            
            .extract_matrix_from_Hessian(i);
            
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    derivative_minor_Hugoniot_matrix(row, col)(i, j) = ;
                }
            }
        }
    }
    */

//    int eq_number = n - 1; // = 2.
/* This below was valid for n = 3.
    int eq_number = 1; // = 2.
    int i = 0;



    for (int ii = 0; ii < n; ii++){
        if (ii != eq_number){
//            H(i) = diff_F(i)*diff_G(eq_number) - diff_F(eq_number)*diff_G(i);

//            for (int j = 0; j < n; j++) nablaH(j, i) = JF(i, j)*diff_G(eq_number) + JG(eq_number, j)*diff_F(i) - 
//                                                       JF(eq_number, j)*diff_G(i) - JG(i, j)*diff_F(eq_number);

            H(i) = diff_F(ii)*diff_G(eq_number) - diff_F(eq_number)*diff_G(ii);

            for (int j = 0; j < n; j++) nablaH(j, i) = JF(ii, j)*diff_G(eq_number) + JG(eq_number, j)*diff_F(ii) - 
                                                       JF(eq_number, j)*diff_G(ii) - JG(ii, j)*diff_F(eq_number);


            i++;
        }
        
    }
    */

    return;
}

RealVector HugoniotContinuation_nD_nm1D::orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane){
    //return vector_product(matrix_column(hyperplane, 0), matrix_column(hyperplane, 1));
        
    int nm1 = hyperplane.rows(); // Number of variables: n - 1

    RealVector Hugoniot_direction(nm1);

    DoubleMatrix A(hyperplane);
    A.resize(nm1, nm1);

    for (int i = 0; i < nm1; i++) A(i, nm1 - 1) = previous_direction(i);

    A.Gram_Schmidt_orthogonalization_by_columns();

    for (int i = 0; i < nm1; i++) Hugoniot_direction(i) = A(i, nm1 - 1);

    return Hugoniot_direction;
}

HugoniotContinuation_nD_nm1D::HugoniotContinuation_nD_nm1D(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb) : HugoniotContinuation(ff, gg, bb){
    there_is_a_bifurcation_space = false;
}

HugoniotContinuation_nD_nm1D::~HugoniotContinuation_nD_nm1D(){
}


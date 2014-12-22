#include "HugoniotContinuation_nDnD.h"

// TODO: Move this method upwards, into SubPhysics.
//
//void HugoniotContinuation_nDnD::jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH){
void HugoniotContinuation_nDnD::jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                                             const RealVector &G, const DoubleMatrix &JG,
                                             const RealVector &C, const DoubleMatrix &JC,
                                             RealVector &H, DoubleMatrix &nablaH){

    // The number of columns of the Jacobians is the dimension of the space.
    int n     = JF.cols();

    H.resize(n - 1);
    nablaH.resize(n, n - 1);

    // [F] & [G]
    //
    RealVector diff_F = F - ref.F;
    RealVector diff_G = G - ref.G;

    // TODO: Find out why these lines below only work for eq_number = 2.
    //       URGENT!!! DUE YESTERDAY!!!
    //
//    int eq_number = n - 1; // = 2.
    int eq_number = 0; // Equation number MUST be between 0 and (n - 1), or else a segfault will happen.
//    int i = 0;

    int neqm1 = F.size() - 1;
    std::vector<int> index(neqm1);
    for (int i = 0; i < neqm1; i++){
        if (i < eq_number) index[i] = i;
        else               index[i] = i + 1;
    }

//    for (int ii = 0; ii < nm1; ii++){
//        if (ii != eq_number){
//            H(i) = diff_F(i)*diff_G(eq_number) - diff_F(eq_number)*diff_G(i);

//            for (int j = 0; j < n; j++) nablaH(j, i) = JF(i, j)*diff_G(eq_number) + JG(eq_number, j)*diff_F(i) - 
//                                                       JF(eq_number, j)*diff_G(i) - JG(i, j)*diff_F(eq_number);

//            i++;
//        }
//        
//    }

    double diff_F_eq = diff_F(eq_number);
    double diff_G_eq = diff_G(eq_number);

    for (int ii = 0; ii < index.size(); ii++){
        int i = index[ii];

            H(ii) = diff_F(i)*diff_G_eq - diff_F_eq*diff_G(i);

            for (int j = 0; j < n; j++) nablaH(j, ii) = JF(i, j)*diff_G_eq + JG(eq_number, j)*diff_F(i) - 
                                                       JF(eq_number, j)*diff_G(i) - JG(i, j)*diff_F_eq;
    }

    for (int ii = 0; ii < C.size(); ii++){
            H(ii + index.size()) = C(ii);

            for (int j = 0; j < n; j++) nablaH(j, ii + index.size()) = JC(ii, j);
    }

//    double norm_diff_G_squared = diff_G*diff_G;
//    double inner_prod_diff_F_diff_G = diff_F*diff_G;

//    for (int k = 0; k < m; k++){
//        H(k) = diff_F(k)*norm_diff_G_squared - diff_G(k)*inner_prod_diff_F_diff_G;

//        for (int j = 0; j < n; j++){
//            double inner_prod_diff_G_column_JG = 0.0;
//            for (int i = 0; i < m; i++) inner_prod_diff_G_column_JG += diff_G(i)*JG(, );

//            nablaH(j, k) = JF(j, k)*norm_diff_G_squared + 2.0*diff_F(k)*
//        }
//    }

    return;
}

RealVector HugoniotContinuation_nDnD::orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane){
    int n = hyperplane.rows();

    RealVector Hugoniot_direction(n);

    DoubleMatrix A(hyperplane);
    A.resize(n, n);

    for (int i = 0; i < n; i++) A(i, n - 1) = previous_direction(i);

    A.Gram_Schmidt_orthogonalization_by_columns();

    for (int i = 0; i < n; i++) Hugoniot_direction(i) = A(i, n - 1);

    return Hugoniot_direction;
}

// TODO: Move this method upwards, into SubPhysics.
//
//int HugoniotContinuation_nDnD::fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &hyperplane, RealVector &Hugoniot_direction){
//    Hugoniot_direction = orthogonalize(previous_direction, hyperplane); // Use Gram-Schmidt BY COLUMNS. A = [nablaH | previous]. Don't normalize (it's done below).
//    
//    double norm_Hugoniot_direction = norm(Hugoniot_direction);
//    
//    // TODO: The 1e-4 below should be variable and should come from outside somehow.
//    //
//    if (norm_Hugoniot_direction < 1e-4) return HUGONIOTCONTINUATION_DIRECTION_ERROR;
//    else {
//        double inv = 1.0/norm_Hugoniot_direction;
//        Hugoniot_direction = Hugoniot_direction*inv;
//    
//        if (Hugoniot_direction*previous_direction < 0.0) Hugoniot_direction = -Hugoniot_direction;
//
//        return HUGONIOTCONTINUATION_DIRECTION_OK;
//    }
//}

HugoniotContinuation_nDnD::HugoniotContinuation_nDnD(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb) : HugoniotContinuation(ff, gg, bb) {
    there_is_a_bifurcation_space = false;
}

HugoniotContinuation_nDnD::~HugoniotContinuation_nDnD(){
}


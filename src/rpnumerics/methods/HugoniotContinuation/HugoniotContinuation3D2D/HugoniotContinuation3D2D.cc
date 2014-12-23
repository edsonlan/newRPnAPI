#include "HugoniotContinuation3D2D.h"

//void HugoniotContinuation3D2D::jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH){
void HugoniotContinuation3D2D::jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                                            const RealVector &G, const DoubleMatrix &JG,
                                            const RealVector &C, const DoubleMatrix &JC,
                                            RealVector &H, DoubleMatrix &nablaH){
//    int n = p.size();

    // The number of columns of the Jacobians is the dimension of the space.    
    int n = JF.cols();

    H.resize(n - 1); // Was: n
    nablaH.resize(n, n - 1);

//    RealVector F(n), G(n);
//    DoubleMatrix JF(n, n), JG(n, n);

//    f->fill_with_jet(n, p.components(), 1, F.components(), JF.data(), 0);
//    g->fill_with_jet(n, p.components(), 1, G.components(), JG.data(), 0);

    // [F] & [G]
    //
    RealVector diff_F = F - ref.F;
    RealVector diff_G = G - ref.G;

//    int eq_number = n - 1; // = 2.
    int eq_number = 1; // = 2.
    int i = 0;

//    for (int i = 0; i < n - 1; i++){
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

    return;
}

RealVector HugoniotContinuation3D2D::orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane){
    return vector_product(matrix_column(hyperplane, 0), matrix_column(hyperplane, 1));
}

//int HugoniotContinuation3D2D::fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &hyperplane, RealVector &Hugoniot_direction){
//
//    //Hugoniot_direction = vector_product(matrix_column(hyperplane, 0), matrix_column(hyperplane, 1));
//    Hugoniot_direction = orthogonalize(previous_direction, hyperplane);
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

HugoniotContinuation3D2D::HugoniotContinuation3D2D(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb) : HugoniotContinuation(ff, gg, bb){
    there_is_a_bifurcation_space = false;
}

HugoniotContinuation3D2D::~HugoniotContinuation3D2D(){
}


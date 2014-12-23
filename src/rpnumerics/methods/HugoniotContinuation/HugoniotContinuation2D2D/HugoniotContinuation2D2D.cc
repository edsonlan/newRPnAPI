#include "HugoniotContinuation2D2D.h"

// TODO: Move this method upwards, into SubPhysics.
//
//void HugoniotContinuation2D2D::jet_Hugoniot(const RealVector &p, RealVector &H, DoubleMatrix &nablaH){
void HugoniotContinuation2D2D::jet_Hugoniot(const RealVector &F, const DoubleMatrix &JF, 
                                            const RealVector &G, const DoubleMatrix &JG,
                                            const RealVector &C, const DoubleMatrix &JC,
                                            RealVector &H, DoubleMatrix &nablaH){
    int n = 2;

//    RealVector F(n), G(n);
//    DoubleMatrix JF(n, n), JG(n, n);

//    f->fill_with_jet(n, p.components(), 1, F.components(), JF.data(), 0);
//    g->fill_with_jet(n, p.components(), 1, G.components(), JG.data(), 0);

    // [F] & [G]
    //
    RealVector diff_F = F - ref.F;
    RealVector diff_G = G - ref.G;

    H.resize(1);
    H(0) = diff_F(0)*diff_G(1) - diff_F(1)*diff_G(0);

    // nablaH = 2 x 1
    nablaH.resize(2, 1);
    for (int j = 0; j < 2; j++) nablaH(j, 0) = JF(0, j)*diff_G(1) + JG(1, j)*diff_F(0) - JF(1, j)*diff_G(0) - JG(0, j)*diff_F(1);

    return;
}

RealVector HugoniotContinuation2D2D::orthogonalize(const RealVector &previous_direction, const DoubleMatrix &hyperplane){
    RealVector Hugoniot_direction(2);
    
    Hugoniot_direction(0) =  hyperplane(1, 0);
    Hugoniot_direction(1) = -hyperplane(0, 0);
    
    return Hugoniot_direction;
}

// TODO: Move this method upwards, into SubPhysics.
//
//int HugoniotContinuation2D2D::fill_Hugoniot_direction(const RealVector &previous_direction, const DoubleMatrix &subspace, RealVector &Hugoniot_direction){
//
//    Hugoniot_direction.resize(2);
//    Hugoniot_direction(0) =  subspace(1, 0);
//    Hugoniot_direction(1) = -subspace(0, 0);
//    
//    // TODO: The lines above should become orthogonalize(). fill_Hugoniot_direction
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

HugoniotContinuation2D2D::HugoniotContinuation2D2D(const FluxFunction *ff, const AccumulationFunction *gg, const Boundary *bb) : HugoniotContinuation_nDnD(ff, gg, bb){
}

HugoniotContinuation2D2D::~HugoniotContinuation2D2D(){
}


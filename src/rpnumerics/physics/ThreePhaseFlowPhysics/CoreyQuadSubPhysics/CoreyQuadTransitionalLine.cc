#include "CoreyQuadTransitionalLine.h"
#include "CoreyQuadSubPhysics.h"

CoreyQuadTransitionalLine::CoreyQuadTransitionalLine(CoreyQuadSubPhysics *c): BifurcationCurve(), corey_(c){
}

CoreyQuadTransitionalLine::~CoreyQuadTransitionalLine(){
}

//double CoreyQuadTransitionalLine::evaluate_point(int type, const RealVector &p){
//    double sw = p(0);
//    double so = p(1);

//    double muw = corey_->muw()->value();
//    double muo = corey_->muo()->value();
//    double mug = corey_->mug()->value();

//    if (type == COREYQUADTRANSITIONALLINE_O_B){
//        return muw*(1.0 - so) - (muw + mug)*sw;
//    }
//    else if (type == COREYQUADTRANSITIONALLINE_W_E){
//        return (mug + muo)*so - muo*(1.0 - sw);
//    }
//    else if (type == COREYQUADTRANSITIONALLINE_G_D){
//        return muw*so - muo*sw;
//    }
//}

double CoreyQuadTransitionalLine::WE_line(void *obj, const RealVector &p){
    CoreyQuadTransitionalLine *cqtl = (CoreyQuadTransitionalLine*)obj;

    double sw = p(0);
    double so = p(1);

    double muw = cqtl->corey_->muw()->value();
    double muo = cqtl->corey_->muo()->value();
    double mug = cqtl->corey_->mug()->value();

    return (mug + muo)*so - muo*(1.0 - sw);
}

double CoreyQuadTransitionalLine::OB_line(void *obj, const RealVector &p){
    CoreyQuadTransitionalLine *cqtl = (CoreyQuadTransitionalLine*)obj;

    double sw = p(0);
    double so = p(1);

    double muw = cqtl->corey_->muw()->value();
    double muo = cqtl->corey_->muo()->value();
    double mug = cqtl->corey_->mug()->value();

    return muw*(1.0 - so) - (muw + mug)*sw;
}

double CoreyQuadTransitionalLine::GD_line(void *obj, const RealVector &p){
    CoreyQuadTransitionalLine *cqtl = (CoreyQuadTransitionalLine*)obj;

    double sw = p(0);
    double so = p(1);

    double muw = cqtl->corey_->muw()->value();
    double muo = cqtl->corey_->muo()->value();
    double mug = cqtl->corey_->mug()->value();

    return muw*so - muo*sw;
}


void CoreyQuadTransitionalLine::curve(int type, std::vector<RealVector> &c){
    c.clear();

    if (type == COREYQUADTRANSITIONALLINE_O_B){
        c.push_back(corey_->O());
        c.push_back(corey_->B());
    }
    else if (type == COREYQUADTRANSITIONALLINE_W_E){
        c.push_back(corey_->W());
        c.push_back(corey_->E());
    }
    else if (type == COREYQUADTRANSITIONALLINE_G_D){
        c.push_back(corey_->G());
        c.push_back(corey_->D());
    }

    return;
}


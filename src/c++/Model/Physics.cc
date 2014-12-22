#include "Physics.h"

Physics::~Physics(){
    // A given subphysics could be instantiated depending on a previously
    // instantiated subphysics. Therefore, they are destroyed in
    // reverse order.
    //
    for (int i = subphysics_.size() - 1; i >= 0; i--) delete subphysics_[i];

    // A given auxiliary function could be instantiated depending on a previously
    // instantiated auxiliary function. Therefore, they are destroyed in
    // reverse order.
    //
    // Subphysics are destroyed first because they are expected to depend
    // on the auxiliary functions.
    // 
    for (int i = auxiliaryfunction_.size(); i >=0 ; i--) delete auxiliaryfunction_[i];
}


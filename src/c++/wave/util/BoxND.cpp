#include "BoxND.h"

BoxND::BoxND(){
}

BoxND::BoxND(const PointND &ppmin, const PointND &ppmax){
//BoxND::BoxND(PointND &ppmin, PointND &ppmax){
    int n = ppmin.size();

    pmin = PointND(n);//ppmin;
    pmax = PointND(n);//ppmax;

    for (int i = 0; i < n; i++){
        pmin.component(i) = min(ppmin(i), ppmax(i));
        pmax.component(i) = max(ppmin(i), ppmax(i));
    }
}

BoxND::BoxND(const BoxND &b){
    pmin = b.pmin;
    pmax = b.pmax;
}

BoxND::BoxND(const BoxND *b){
    pmin = b->pmin;
    pmax = b->pmax;
}

BoxND::~BoxND(){
}

// Box-box collision test, very similar to circle-circle collision test.
// Verify that the distance between the centers of the boxes is less
// than the sum of the "radii" for each dimension.
//
bool BoxND::intersect(const BoxND &b) const {
    int n = pmin.size();

    PointND bpmin = b.pmin, bpmax = b.pmax;

    double this_radius[n], b_radius[n];

    for (int i = 0; i < n; i++){
        this_radius[i] = (pmax.component(i)  - pmin.component(i))/2.0;
        b_radius[i]    = (bpmax.component(i) - bpmin.component(i))/2.0;
    }

    double this_center[n], b_center[n];

    for (int i = 0; i < n; i++){
        this_center[i] = pmin.component(i) + this_radius[i];
        b_center[i]    = bpmin.component(i) + b_radius[i];
    }

    bool collision = true; int pos = 0;

    while (collision && pos < n){
        if (abs(this_center[pos] - b_center[pos]) > (this_radius[pos] + b_radius[pos])) collision = false;
        pos++;
    }

    return collision;
}


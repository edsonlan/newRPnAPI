/* 
 * File:   ViewList.cpp
 * Author: edsonlan
 * 
 * Created on September 29, 2014, 3:53 PM
 */

#include "ViewList.h"

ViewList::ViewList() : curveViewList_(new vector<RpnSolution *>()) {
}

void ViewList::add(RpnSolution * curveView) {

    curveViewList_->push_back(curveView);

}

RpnSolution * ViewList::get(int n) {


    return curveViewList_->at(n);



}

int ViewList::getSize() {

    return curveViewList_->size();

}

void ViewList::clean() {


    for (int i = 0; i < curveViewList_->size(); i++) {


        delete curveViewList_->at(i);

    }


    curveViewList_->clear();



}

void ViewList::remove(RpnSolution * curveView) {

    for (std::vector<RpnSolution *>::iterator it = curveViewList_->begin(); it != curveViewList_->end(); ++it) {

        if (*it == curveView) {
            curveViewList_->erase(it);
            delete curveView;
        }
    }




}

ViewList::ViewList(const ViewList& orig) {
}

ViewList::~ViewList() {

    for (int i = 0; i < curveViewList_->size(); i++) {

        delete get(i);

    }


    delete curveViewList_;
}


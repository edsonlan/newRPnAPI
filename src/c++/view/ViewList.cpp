/* 
 * File:   ViewList.cpp
 * Author: edsonlan
 * 
 * Created on September 29, 2014, 3:53 PM
 */

#include "ViewList.h"

ViewList::ViewList() : curveViewList_(new vector<RPnResult *>()) {
}

void ViewList::add(RPnResult * curveView) {

    curveViewList_->push_back(curveView);

}

RPnResult * ViewList::get(int n) {


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

void ViewList::remove(RPnResult * curveView) {

    for (std::vector<RPnResult *>::iterator it = curveViewList_->begin(); it != curveViewList_->end(); ++it) {

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


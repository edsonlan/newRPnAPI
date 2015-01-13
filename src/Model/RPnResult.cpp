/* 
 * File:   Result.cpp
 * Author: edsonlan
 * 
 * Created on December 16, 2014, 1:09 PM
 */

#include "RPnResult.h"
#include "RPnMethod.h"


RPnResult::RPnResult( RPnMethod * method, Data * coords) : method_(method), data_(coords) {
}

RPnResult::RPnResult(const RPnResult& orig) : method_(orig.method_) {
}

void RPnResult::recalc() {

    method_->recalc(data_);

}

RPnResult::~RPnResult() {

    delete data_;
    delete method_;


}

const Data  & RPnResult::getData() const {
    return *data_;
}

const RPnMethod & RPnResult::getMethod() const {

    return *method_;
}
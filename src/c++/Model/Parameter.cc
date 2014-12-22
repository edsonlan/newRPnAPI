#include "Parameter.h"

void Parameter::copy(const Parameter *original){
    name_  = original->name_;
    min_   = original->min_;
    max_   = original->max_;
    value_ = original->value_;

    return;
}

Parameter::Parameter() : name_(std::string("")), 
                         min_(0.0),
                         max_(0.0),
                         value_(0.0){
}

Parameter::Parameter(const std::string &n, double v) : name_(n),
                                                       min_(v),
                                                       max_(v),
                                                       value_(v) {
}

Parameter::Parameter(const std::string &n, double mn, double mx, double v) : name_(n),
                                                                             min_(mn),
                                                                             max_(mx),
                                                                             value_(v) {
}

Parameter::Parameter(const Parameter &original){
    copy(&original);
}

Parameter::Parameter(const Parameter *original){
    copy(original);
}

Parameter::~Parameter(){
}

Parameter Parameter::operator=(const Parameter &original){
    if (this != &original) copy(&original);

    return *this;
}


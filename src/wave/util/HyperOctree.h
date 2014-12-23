#ifndef _HYPEROCTREE_
#define _HYPEROCTREE_

#include <sstream>
#include <string>
#include <vector>
#include <set>
#include "BoxND.h"

// Class HyperOctree.
//
// This class implements a hyperoctree, of arbitrary dimension.

template <typename T>
class HyperOctree {
    private:
    protected:
        BoxND box_;        

        unsigned long level_, maxlevel_, maxsize_;

        HyperOctree<T> *father_;
        std::vector<HyperOctree<T>*> sons_;

        std::set<T*> data_;

        void subdivide_box(BoxND &b, std::vector<BoxND> &vb);

        // Recover the data stored on the current node's sons.
        // Returns true if the list has up to maxsize_ elements,
        // false otherwise. If true, the sons can be deleted and
        // the current node can become a leaf.
        //
        // By returning a boolean, it is not necessary to recover
        // all the data stored in the sons, thus shortening the process.
        //
        bool recover_list(std::set<T*> &list);

        void copy(const HyperOctree<T> *orig);

        // The engine for the structure.
        //
        std::string structure(std::string &path, std::string &s) const ;
    public:
        HyperOctree(const BoxND &b, unsigned long maxlevel = 3, unsigned long maxsize = 100);
        HyperOctree(const HyperOctree<T> *orig);
        HyperOctree(const HyperOctree<T> &orig);

        ~HyperOctree();

        void add(T* obj);
        void remove(T* obj);
        void clear();

        void within_box(const BoxND &b, std::set<T*> &list) const;
        void within_box(const BoxND &b, std::vector<T*> &list) const;

        unsigned long number_of_objects(void) const;

        std::string structure() const;

        void boxes(std::vector<BoxND> &boxes_of_nodes, std::vector<BoxND> &boxes_of_leaves) const;

        void leaves_intersected_by_object(T *object, std::vector<BoxND> &boxes) const;

        void within_leaves_intersected_by_object(T *object, std::set<T*> &list) const;
        void within_leaves_intersected_by_object(T *object, std::vector<T*> &list) const;
};

// Subdivide a box unsigned longo 2^space_dimension sub-boxes.
//
template <typename T>
void HyperOctree<T>::subdivide_box(BoxND &b, std::vector<BoxND> &vb){
    vb.clear();

    unsigned long n = b.pmin.size();

    PointND ppmin(n), ppmax(n);
    double mid = .5*(b.pmin.component(0) + b.pmax.component(0));

    if (n == 1){
        ppmin.component(0) = b.pmin.component(0);
        ppmax.component(0) = mid;
        vb.push_back(BoxND(ppmin, ppmax));

        ppmin.component(0) = mid;
        ppmax.component(0) = b.pmax.component(0);
        vb.push_back(BoxND(ppmin, ppmax));
    }
    else {
        PointND tpmin(n - 1), tpmax(n - 1);
        for (unsigned long i = 1; i < n; i++){
            tpmin.component(i - 1) = b.pmin.component(i);
            tpmax.component(i - 1) = b.pmax.component(i);
        }

        BoxND tb(tpmin, tpmax);
        std::vector<BoxND> tvb;
        subdivide_box(tb, tvb);

        for (unsigned long i = 0; i < tvb.size(); i++){
            for (unsigned long j = 1; j < n; j++){
                ppmin.component(j) = tvb[i].pmin.component(j - 1);
                ppmax.component(j) = tvb[i].pmax.component(j - 1);
            }

            ppmin.component(0) = b.pmin.component(0);
            ppmax.component(0) = mid;
            vb.push_back(BoxND(ppmin, ppmax));

            ppmin.component(0) = mid;
            ppmax.component(0) = b.pmax.component(0);
            vb.push_back(BoxND(ppmin, ppmax));
        }
    }

    return;
}

// Recover the data stored on the current node's sons.
// Returns true if the list has up to maxsize_ elements,
// false otherwise. If true, the sons can be deleted and
// the current node can become a leaf.
//
// By returning a boolean, it is not necessary to recover
// all the data stored in the sons, thus shortening the process.
//
template <typename T>
bool HyperOctree<T>::recover_list(std::set<T*> &list){
    // If leaf
    if (sons_.size() == 0){
        if (data_.size() + list.size() > maxsize_) return false;
        else {
            typename std::set<T*>::iterator it;
            for (it = data_.begin(); it != data_.end(); it++) list.insert(*it);

            return true;
        }
    }
    // If node
    else {
        bool valid = true;
        unsigned long i = 0;

        while (valid && i < sons_.size()){
            valid = sons_[i]->recover_list(list);
            i++;
        }

        return valid;
    }
}

template <typename T>
HyperOctree<T>::HyperOctree(const BoxND &b, unsigned long maxlevel, unsigned long maxsize) :  
                            box_(b), maxlevel_(maxlevel), maxsize_(maxsize), level_(0), father_(0){
}

template <typename T>
void HyperOctree<T>::copy(const HyperOctree<T> *orig){
    if (orig->sons_.size() == 0){
        typename std::set<T*>::iterator it;
        for (it = orig->data_.begin(); it != orig->data_.end(); it++) data_.insert(*it);
    }
    else {
        sons_.resize(orig->sons_.size());
        for (unsigned long i = 0; i < orig->sons_.size(); i++){
             HyperOctree *temp = new HyperOctree(orig->sons_[i]);
             temp->father_ = this;
             temp->level_ = level_ + 1;
             sons_[i] = temp;
        }
    }

    return;
}

template <typename T>
HyperOctree<T>::HyperOctree(const HyperOctree<T> *orig){
    box_      = orig->box_;
    maxlevel_ = orig->maxlevel_;
    maxsize_  = orig->maxsize_;
    level_    = 0;
    father_   = 0;

    copy(orig);
}

template <typename T>
HyperOctree<T>::HyperOctree(const HyperOctree<T> &orig){
    box_      = orig.box_;
    maxlevel_ = orig.maxlevel_;
    maxsize_  = orig.maxsize_;
    level_    = 0;
    father_   = 0;

    copy(&orig);
}

template <typename T>
HyperOctree<T>::~HyperOctree(){
    clear();
}

template <typename T>
void HyperOctree<T>::add(T* obj){
    if (obj->intersect(box_)){
        // If leaf
        if (sons_.size() == 0){
            data_.insert(obj);     

            // Subdivide if need be.
            //
            if (data_.size() > maxsize_ && level_ < maxlevel_){
                std::vector<BoxND> vb;
                subdivide_box(box_, vb);

                // Create the sons... 
                for (unsigned long i = 0; i < vb.size(); i++){
                    HyperOctree *temp = new HyperOctree(vb[i], maxlevel_, maxsize_);
                    temp->level_ = level_ + 1;
                    temp->father_ = this;
    
                    // ...redistribute the data...
                    typename std::set<T*>::iterator it;
                    for (it = data_.begin(); it != data_.end(); it++) temp->add(*it);

                    sons_.push_back(temp);
                }

                // ...and become a node.
                data_.clear();
            }
        }
        // If node
        else {
            for (unsigned long i = 0; i < sons_.size(); i++) sons_[i]->add(obj);
        }
    }

    return;
}

template <typename T>
void HyperOctree<T>::remove(T* obj){
    if (obj->intersect(box_)){
        // If leaf
        if (sons_.size() == 0){
            data_.erase(data_.find(obj));
        }
        // If node
        else {
            for (unsigned long i = 0; i < sons_.size(); i++) sons_[i]->remove(obj);

            if (recover_list(data_)){
                for (unsigned long i = 0; i < sons_.size(); i++) delete sons_[i];
                sons_.clear();
            }
            else data_.clear();
        }
    }

    return;
}

template <typename T>
void HyperOctree<T>::clear(void){
    data_.clear();

    for (unsigned long i = 0; i < sons_.size(); i++) delete sons_[i];
    sons_.clear();
}

template <typename T>
void HyperOctree<T>::within_box(const BoxND &b, std::set<T*> &list) const {
    if (father_ == 0) list.clear();

    if (box_.intersect(b)){
        // If leaf:
        if (sons_.size() == 0){
            typename std::set<T*>::iterator it;
            for (it = data_.begin(); it != data_.end(); it++) if ((*it)->intersect(b)) list.insert(*it);
        }
        // If node:
        else {
            for (unsigned long i = 0; i < sons_.size(); i++) sons_[i]->within_box(b, list);
        }
    }

    return;
}

template <typename T>
void HyperOctree<T>::within_box(const BoxND &b, std::vector<T*> &list) const {
    std::set<T*> list_set;
    within_box(b, list_set);

    list.resize(list_set.size());

    typename std::set<T*>::iterator it; unsigned long i;
    for (it = list_set.begin(), i = 0; it != list_set.end(); it++, i++) list[i] = *it;

    return;
}

template <typename T>
unsigned long HyperOctree<T>::number_of_objects(void) const {
    unsigned long n = 0;

    // If leaf:
    if (sons_.size() == 0) n = data_.size();
    // If node:
    else {
        for (unsigned long i = 0; i < sons_.size(); i++) n += sons_[i]->number_of_objects();
    }

    return n;
}

// Protected, engine.
//
template <typename T>
std::string HyperOctree<T>::structure(std::string &path, std::string &s) const {
    // Update the path
    std::string this_path(path);
    if (father_ != 0){
        unsigned long i = 0; bool found = false;
        while (this != father_->sons_[i]) i++;

        std::stringstream sspath;
        sspath << i;

        this_path += "." + sspath.str();
    }

    std::string str(s);
    
    if (father_ == 0) str += "Root. ";
    else str += "--";

    str += (sons_.size() != 0) ? "Node." : "Leaf.";
    str += " Path: " + this_path + ". Box: ";

    std::stringstream ss;
    ss << box_.pmin << "-" << box_.pmax;
    str += ss.str();

    if (sons_.size() != 0) str += "\n "; // str += "\n>";

    std::string next(s);
    if (father_ != 0 && this == father_->sons_[father_->sons_.size() - 1]) next.replace(next.rfind("|"), 1, " ");

    // If leaf:
    if (sons_.size() == 0){
        typename std::set<T*>::iterator it;
        for (it = data_.begin(); it != data_.end(); it++){
            str += "\n "; // str += "\n&";
            str += next + "  |-";
            std::stringstream temp;
            temp << "Pointer: " << *it;
            str += temp.str();
        }
        str += "\n " + next;
    }
    // If node
    else {
        //str += "\n@";
        next += "  |";
        for (unsigned long i = 0; i < sons_.size(); i++) str += sons_[i]->structure(this_path, next);
    }

    if (father_ != 0 && this != father_->sons_[father_->sons_.size() - 1]) str += "\n "; //str += "\n" + next + "\n*";

    return str;    
}

// Public unsigned longerface.
//
template <typename T>
std::string HyperOctree<T>::structure() const {
    std::string path("0");
    std::string s("");

    return structure(path, s);
}

template <typename T>
void HyperOctree<T>::boxes(std::vector<BoxND> &boxes_of_nodes, std::vector<BoxND> &boxes_of_leaves) const {
    if (father_ == 0){
        boxes_of_nodes.clear();
        boxes_of_leaves.clear();
    }

    if (sons_.size() == 0) boxes_of_leaves.push_back(box_);
    else {
        boxes_of_nodes.push_back(box_);
        for (unsigned long i = 0; i < sons_.size(); i++) sons_[i]->boxes(boxes_of_nodes, boxes_of_leaves);
    }

    return;
}

template <typename T>
void HyperOctree<T>::leaves_intersected_by_object(T *object, std::vector<BoxND> &boxes) const {
    if (father_ == 0) boxes.clear();

    if (object->intersect(box_)){
        if (sons_.size() == 0) boxes.push_back(box_);
        else {
            for (unsigned long i = 0; i < sons_.size(); i++) sons_[i]->leaves_intersected_by_object(object, boxes);
        }
    }

    return;
}

template <typename T>
void HyperOctree<T>::within_leaves_intersected_by_object(T *object, std::set<T*> &list) const {
    if (father_ == 0) list.clear();

    if (object->intersect(box_)){
        // If leaf:
        if (sons_.size() == 0){
            typename std::set<T*>::iterator it;
            for (it = data_.begin(); it != data_.end(); it++) list.insert(*it);
        }
        // If node:
        else {
            for (unsigned long i = 0; i < sons_.size(); i++) sons_[i]->within_leaves_intersected_by_object(object, list);
        }
    }

    return;
}

template <typename T>
void HyperOctree<T>::within_leaves_intersected_by_object(T *object, std::vector<T*> &list) const {
    std::set<T*> list_set;
    within_leaves_intersected_by_object(object, list_set);

    list.resize(list_set.size());

    typename std::set<T*>::iterator it; unsigned long i;
    for (it = list_set.begin(), i = 0; it != list_set.end(); it++, i++) list[i] = *it;

    return;
}

#endif // _HYPEROCTREE_


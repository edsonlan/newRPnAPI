#ifndef _MULTIARRAY_
#define _MULTIARRAY_

#include <vector>
#include <iostream>
#include <cstdarg>

// Class MultiArray implements an array of doubles of arbitrary dimension_.
// Elements within MultiArray are indexed by a vector, also of size dimension_.
// Internally, a std::vector<double> is used.
//
// Every dimension of a MultiArray is constrained
// to a range. Suppose, without losing generality, that a 3D space is dealt with
// and let the ranges be Nx, Ny, Nz. Thus, an element [i j k] is such that
//
//     0 <= i < Nx, 0 <= j < Ny, 0 <= k < Nz.
//
// (The outmost index is i, the innermost is k.) Suppose that
// range_ = [2 3 4]. Then the elements are laid out in memory
// as:
//
//                   k
//             0   1   2   3             
//
//          |  0   1   2   3 |   0
// i = 0 => |  4   5   6   7 |   1  j
//          |  8   9  10  11 |   2
//
// and
//                   k
//             0   1   2   3
//
//          | 12  13  14  15 |  0
// i = 1 => | 16  17  18  19 |  1  j
//          | 20  21  22  23 |  2
//
// and multiindex [i j k] points to the storage position
// 
//     p = i*(Nj*Nk) + j*Nk + k*1.
//
// This direct mapping can be conveniently written as the
// inner product of two vectors.
//
// The reverse mapping, that is, finding the multiindex given a location 
// in the array, is implemented by withdrawing
// the contribution of each dimension from the location
// (whenever the symbol "/" is used, assume integer division):
//
//     i = p/(Nj*Nk),
//     j = (p - i*(Nj*Nk))/Nk,
//     k = (p - i*(Nj*Nk) - j*Nk)/1.
//
// Observe that in both cases it is convenient to have
// the factors Nj*Nk, etc., computed beforehand and readily
// available.
//
template <typename T>
class MultiArray {
    private:
    protected:
        // The dimension of the MultiArray, dimension_ = range_.size();
        //
        int dimension_;

        // Where the data are stored and the size.
        //
        T *storage_;

        long unsigned int size_;

        // The range of each dimension of the MultiArray.
        //
        std::vector<long unsigned int> range_;

        // Auxiliary vector used for mapping a multiindex into a location
        // and viceversa.
        //
        std::vector<long unsigned int> mult;

        // Copy a MultiArray (used by the copy-constructors).
        //
        void copy(const MultiArray<T> *orig);
    public:
        MultiArray();
        MultiArray(const std::vector<unsigned long int> &orig_range_);
        MultiArray(const MultiArray<T> *orig);
        MultiArray(const MultiArray<T> &orig);

        virtual ~MultiArray();

        MultiArray<T> operator=(const MultiArray<T> &original);

        // Resize a MultiArray. The new dimension is orig_range_.size().
        // This method is useful because std::vector<MultiArray> and
        // other containers could thus be used.
        //
        virtual void resize(const std::vector<unsigned long int> &orig_range_);

        // The dimension of the MultiArray.
        //
        virtual int dimension() const {
            return dimension_;
        }

        // The size of the MultiArray, that is,
        //
        //     range_[0]*range_[1]*...*range_[dimension_ - 1].
        //
        virtual unsigned long int size() const {
            return size_;
        }

        virtual std::vector<long unsigned int> range() const {
            return range_;
        }

        // Convert a location into a multiindex.
        //
        virtual std::vector<unsigned long int> multiindex(unsigned long int location) const;

        // Convert a multiindex into a location (classical).
        //
        virtual unsigned long int location(const std::vector<unsigned long int> &multiindex) const;

        // Convert a multiindex into a location (variadic).
        //
        virtual unsigned long int location(...) const;

        // Access a position in the array (classical).
        //
        virtual const T& operator()(const std::vector<unsigned long int> &multiindex) const;
        virtual       T& operator()(const std::vector<unsigned long int> &multiindex);

        // Access a position in the array (variadic).
        //
        virtual const T& operator()(...) const;
        virtual       T& operator()(...);

        // Ouput to a stream.
        //
        friend std::ostream & operator<<(std::ostream &os, const MultiArray<T> &ma){
            std::vector<unsigned long int> mi;

            for (int i = 0; i < ma.size_; i++){
                mi = ma.multiindex(i);

                os << "[";
                for (int j = 0; j < mi.size(); j++){
                    os << mi[j];

                    if (j < mi.size() - 1) os << ", ";
                }
                os << "] = " << ma.storage_[i] << std::endl;
            }

            return os;
        }
};

template <typename T>
MultiArray<T>::MultiArray(){
    storage_ = 0;
}

template <typename T>
MultiArray<T>::MultiArray(const std::vector<unsigned long int> &orig_range_){
    storage_ = 0;
    resize(orig_range_);
}

template <typename T>
MultiArray<T>::MultiArray(const MultiArray<T> *orig){
    storage_ = 0;
    copy(orig);
}

template <typename T>
MultiArray<T>::MultiArray(const MultiArray<T> &orig){
    storage_ = 0;
    copy(&orig);
}

template <typename T>
MultiArray<T>::~MultiArray(){
    if (storage_ != 0) delete [] storage_;
}

template <typename T>
MultiArray<T> MultiArray<T>::operator=(const MultiArray<T> &original){
    if (this != &original) copy(&original);

    return *this;
}

template <typename T>
void MultiArray<T>::copy(const MultiArray<T> *orig){
    dimension_ = orig->dimension_;
    size_      = orig->size_;
    range_     = orig->range_;
    mult       = orig->mult;

    if (storage_ != 0) delete [] storage_;
    storage_   = new T[size_];
    for (unsigned int i = 0; i < size_; i++) storage_[i] = orig->storage_[i];

    return;
}

template <typename T>
void MultiArray<T>::resize(const std::vector<unsigned long int> &orig_range_){
    range_ = orig_range_;
    dimension_ = range_.size();

    // Prepare mult and n (the total number of elements).
    //
    size_ = 1;
    mult.resize(dimension_);
    for (int i = 0; i < dimension_; i++){
        mult[dimension_ - i - 1] = size_;
        size_ *= range_[dimension_ - i - 1];
    }

    // The total amount of storage needed is the product
    // of all the elements of range_.
    //
    if (storage_ != 0) delete [] storage_;
    storage_ = new T[size_];

    return;
}

template <typename T>
std::vector<unsigned long int> MultiArray<T>::multiindex(unsigned long int location) const {
    std::vector<unsigned long int> mi(dimension_);

    int p = location;

    for (int i = 0; i < dimension_; i++){
        mi[i] = p/mult[i];

        p -= mi[i]*mult[i];
    }

    return mi;
}

template <typename T>
unsigned long int MultiArray<T>::location(const std::vector<unsigned long int> &multiindex) const {
//    int pos = multiindex.back();
//    int prod = 1;

//    for (int i = dimension_ - 1; i > 0; i--){
//        prod *= range_[i];

//        pos += prod*multiindex[i - 1];
//    }

    unsigned long int pos = 0;
    for (int i = 0; i < dimension_; i++){
        pos += mult[i]*multiindex[i];
    }

    return pos;
}

template <typename T>
unsigned long int MultiArray<T>::location(...) const {
    std::vector<unsigned long int> multiindex;

    va_list ap;
    va_start(ap, NULL); // No named arguments were passed to the function, thus, NULL is passed here.

    for(int i = 0; i < dimension_; i++) multiindex.push_back(va_arg(ap, unsigned long int));

    va_end(ap);

    return location(multiindex);
}

template <typename T>
const T& MultiArray<T>::operator()(const std::vector<unsigned long int> &multiindex) const {
    return storage_[location(multiindex)];
}

template <typename T>
T& MultiArray<T>::operator()(const std::vector<unsigned long int> &multiindex){
    return storage_[location(multiindex)];
}

template <typename T>
const T& MultiArray<T>::operator()(...) const {
    std::vector<unsigned long int> multiindex;

    va_list ap;
    va_start(ap, NULL); // No named arguments were passed to the function, thus, NULL is passed here.

    for(int i = 0; i < dimension_; i++) multiindex.push_back(va_arg(ap, unsigned long int));

    va_end(ap);

    return storage_[location(multiindex)];
}

template <typename T>
T& MultiArray<T>::operator()(...){
    std::vector<unsigned long int> multiindex;

    va_list ap;
    va_start(ap, NULL); // No named arguments were passed to the function, thus, NULL is passed here.

    for(int i = 0; i < dimension_; i++) multiindex.push_back(va_arg(ap, unsigned long int));

    va_end(ap);

    return storage_[location(multiindex)];
}

#endif // _MULTIARRAY_


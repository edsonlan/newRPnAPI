#ifndef _QUICKSORT_
#define _QUICKSORT_

#include <vector>
#include <iostream>

template <typename T>
class QuickSort {
    private:
    protected:
        static bool trivial_order(void*, const T&, const T&);
    public:
        static void quicksort(std::vector<T> &v);
        static void quicksort(std::vector<T> &v, void *obj, bool (*f)(void*, const T&, const T&));
};

template <typename T> 
bool QuickSort<T>::trivial_order(void*, const T &a, const T &b){
    return a < b;
}

// Default quicksort. Elements are ordained according to what the less-than operator
// of T defines.
//
template <typename T> 
void QuickSort<T>::quicksort(std::vector<T> &v){
    quicksort(v, 0, &trivial_order);

    return;
}

// Generic quicksort. The user provides a pointer to a function of the form:
//
//     bool f(void *obj, const T &a, const T &b);
//
// that returns true if a < b, false otherwise. The pointer to void obj
// may be used to pass an object, but it is not required.
//
// If f is null, the less-than operator will be used (default).
//
template <typename T> 
void QuickSort<T>::quicksort(std::vector<T> &v, void *obj, bool (*f)(void*, const T&, const T&)){
    // If the comparison function is not provided default to the less-than operator.
    //
    if (f == 0) f = &trivial_order;

    int n = v.size();

    if (n > 1){
        int pivot_index = n/2;
        T pivot = v[pivot_index];

        std::vector<T> smaller, greater;
        for (int i = 0; i < n; i++){
            if (i != pivot_index){
                if ( (*f)(obj, v[i], v[pivot_index]) ) smaller.push_back(v[i]);
                else                                   greater.push_back(v[i]);
            }
        }

        if (smaller.size() > 1) quicksort(smaller, obj, f);
        if (greater.size() > 1) quicksort(greater, obj, f);

        int pos = 0;
        for (int i = 0; i < smaller.size(); i++){
            v[pos] = smaller[i];
            pos++;
        }

        v[pos] = pivot;
        
        for (int i = 0; i < greater.size(); i++){
            pos++;
            v[pos] = greater[i];
        }
    }    

    return;
}

#endif // _QUICKSORT_


/* IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) Space.h
 **/

#ifndef _Space_H
#define	_Space_H

//! Definition of class Space.

/*!
    
TODO:
NOTE :
    
@ingroup rpnumerics
 */

class Space {
private:
    int dim_;
    const char * name_;


public:
    Space(void);
    Space(const char * name, const int dim);
    Space(const Space &);
    ~Space(void);

    const char * name(void) const;
    const int dim(void) const;
    void name(const char * name);
    void dim(const int dim);


};

inline Space::Space() : dim_(2),name_("Space") {
}

inline Space::Space(const char * name, const int dim)
:dim_(dim), name_(name)
 {
}

inline Space::Space(const Space &copy) : dim_(copy.dim()), name_(copy.name()) {
}

inline Space::~Space() {
}

inline const char * Space::name() const {
    return name_;
}

inline const int Space::dim() const {
    return dim_;
}

inline void Space::name(const char * name) {
    name_ = name;
}

inline void Space::dim(const int dim) {
    dim_ = dim;
}

#endif	/* _Space_H */

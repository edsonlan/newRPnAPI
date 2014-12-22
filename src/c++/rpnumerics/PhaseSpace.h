/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) PhaseSpace.h
 **/


//!
/*!
 *
 * TODO:
 * NOTE :
 *
 * @ingroup rpnumerics
 */


#ifndef _PhasSpace_H
#define	_PhasSpace_H

class PhaseSpace {
    
    
public:
    
    PhaseSpace();
    ~PhaseSpace();
    
    //! PhaseSpace name accessor
    const char *getName() const ;
    
    
    //! PhaseSpace name mutator
    void setName(const char *name);
    
    
    //! PhaseSpace dimension accessor
    const int getDim() const ;
    
    
    //! PhaseSpace dimension mutator
    void setDim(int dim) const;
    
    
private :
    
    const   char * name_;
    const  int dim_;
    
};
inline const char * PhaseSpace::getName{ return name_;}

inline const int  PhaseSpace::getDim { return dim_;}

inline void PhaseSpace::setName(const char * name){ name_=name;}

inline void PhaseSpace::setDim(const int dim) {dim_=dim;}

#endif	/* _PhasSpace_H */


#ifndef _PHYSICS_
#define _PHYSICS_

#include <vector>
#include <string>

// In order to instantiate Physics as a plugin:
//
#ifdef __linux__
    #include <dlfcn.h>
#endif

#ifdef __WIN32
    #include "dlfcn.h" // https://code.google.com/p/dlfcn-win32/
#endif

#include "SubPhysics.h"
#include "AuxiliaryFunction.h"

// TODO: Not sure if Rarefaction, and other stuff like that, need be included here. So far
//       it seems to me that not, and that those things should be included in the GUI, be it
//       Java's or FLTK's.
//
//       Morante.

// TODO: Is it necessary to create a class Auxiliary_Function, to be used by the subphysics?
//       For example, Thermodynamic could be such a class. In that way the auxiliary functions
//       could be queried about their parameters, which could be necessary in the GUI.
//       Now, since the Thermodynamic is common for all the Two Phases subphysics, and it
//       makes no sense to have more than one, at least in that model, I believe that it is
//       correct to have this kind of class refered here. However, that doesn't mean that
//       a similar necessity will not arise in the subphysic. Perhaps something along the lines
//       of "Physics Auxiliary Classes" and "Subphysic Auxiliary Classes" is more correct.
//
//       Morante.

// TODO: If a graph is finally used to control how curves are created or whatever, it must be here.
//
//       Morante.

// Class Physics defines a physical model. Several subphysics will be declared herewithin, and pointers
// to them stored in subphysics_ (a vector). Later on, the GUI (or any other superior interface)
// will access this vector to perform computations on the subphysics.
//
// All methods within this class should be virtual, so the derived classes can
// modify them as need be.
//
class Physics {
    private:
        // Users should not add any private members.
    protected:
        // String containing the Physics' name. This string can be used by the GUI to display it.
        //
        // TODO: Since the name is just a descriptor MAYBE it could be better to transform it into a
        //       static member so that the GUI can query the Physics about its name without instantiating it.
        //       Also, maybe it is better to have a member "name" and a member "description," both std::string's.
        //
        //       Morante.
        //
        std::string info_physics_;

        // All Physics models will have one or more SubPhysics (q.v.). SubPhysics will be instantiated
        // by the Physics constructor and live as long as the Physics, but no longer. Therefore, SubPhysics
        // are expected to be destroyed by the Physics' destructor.
        // 
        std::vector<SubPhysics*> subphysics_;

        // Any object with Parameters (q.v.) that are expected to be modified by the user
        // (by the GUI, for instance) is to be coded as a class derived from AuxiliaryFunction (q.v.).
        // Auxiliary functions common to all the SubPhysics will be created by the Physics constructor,
        // and pointers to them will be stored here.
        //
        // The method auxiliaryfunction() below will grant access to this member.
        //
        std::vector<AuxiliaryFunction*> auxiliaryfunction_;
    public:
        // Typically, the subphysics will be instantiated in the physics constructor. Therefore, it
        // is expected that they will be deleted in the physics' destructor.
        //
        // The destructor is virtual because classes will be derived from this one.
        //
        virtual ~Physics();

        // Information about the physic, which can be displayed somewhere.
        //
        virtual const std::string info_physics(){return info_physics_;}

        // To be used by the interface, so it can refer to the subphysics when
        // computing the curves. This method is not const so that the parameters
        // of each subphysics can be modified by the GUI.
        //
        virtual std::vector<SubPhysics*> & subphysics(){return subphysics_;} 

        // To be used by the interface, so it can refer to the auxiliary functions that
        // are common to all the subphysics (such as the Thermodynamics). 
        // This method is not const so that the parameters
        // of each auxiliary function can be modified by the GUI.
        //
        virtual std::vector<AuxiliaryFunction*> & auxiliaryfunction(){return auxiliaryfunction_;} 
};

// To load/unload a Physics-derived class as a plugin. Refer to
//
//     http://www.faqs.org/docs/Linux-mini/C++-dlopen.html
//
typedef Physics* create_physics();
typedef void destroy_physics(Physics*);


/*

class MyPhysics : public Physics {
    STUFF HERE
};

extern "C" Physics* create() {
    return new MyPhysics;
}

extern "C" void destroy(Physics* p) {
    delete p;
}

*/

#endif // _PHYSICS_


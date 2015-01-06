#ifndef _PARAMETER_
#define _PARAMETER_

#include "Subject.h"
#include <string>

// TODO: Decide if minimum and maximum are the correct solution, since there
//       are magnitudes, such as temperature, that do not have an axiomatic
//       range (like, say, saturation), but an arbitrary one. Perhaps
//       in these cases something like 
//
//           double inf = std::numeric_limits<double>::infinity();
//
//       (which is true infinity and can be positive and negative) should be used.
//
//       Morante.

// Class Parameter. Every Auxiliary function or FluxFunction or the like has
// a number of parameters, which have a name, and minimum, maximum and initial 
// (or default) values. Being thus standardized, the interface can query said
// Functions about their parameters and automatically create an appropriate
// user interface to show and modify the parameters interactively.
//
//     Virtual destructor (because classes may be derived from this one).
//
// All methods within this class should be virtual, so the derived classes can
// modify them as need be.
//
class Parameter:public Subject {
    private:
        // Users should not add any private members.
    protected:
        std::string name_;
        double min_, max_, value_;

        virtual void copy(const Parameter *original);
    public:
        Parameter();
        Parameter(const std::string &n, double v);
        Parameter(const std::string &n, double mn, double mx, double v);
        Parameter(const Parameter &original);
        Parameter(const Parameter *original);

        virtual ~Parameter();

        virtual Parameter operator=(const Parameter &original);

        virtual std::string name() const {return name_;}

        virtual double min() const {return min_;}
        virtual void min(double m) {min_ = m; return;}

        virtual double max() const {return max_;}
        virtual void max(double m) {max_ = m; return;}

        virtual double value() const {return value_;}
        virtual void value(double v) {value_ = v; notify(); return;}
};

#endif // _PARAMETERS_


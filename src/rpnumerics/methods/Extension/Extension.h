#ifndef _EXTENSION_
#define _EXTENSION_

#define EXTENSION_OK    0
#define EXTENSION_ERROR 1

#define EXPLICIT_EXTENSION 0
#define IMPLICIT_EXTENSION 1

#include <vector>
#include <string>
#include "RealVector.h"
#include "Curve.h"

class Extension {
    private:
    protected:
        void divide_curve(void *obj, bool (*f)(void*, const RealVector&), const Curve &original, std::vector<Curve> &out);
    public:
        Extension(){
        }

        virtual ~Extension(){
        }

        virtual int extension(const RealVector &p, RealVector &ext_p) = 0;
        virtual std::string name() const = 0;

        // Return EXPLICIT_EXTENSION or IMPLICIT_EXTENSION that the user
        // may know that the curve being extended and its extension are
        // given as a sequence of points or as segments.
        //
        virtual int extension_type() = 0;

        virtual void extension_curve(const Curve &curve, 
                                     void *fdo, bool (*fd)(void*, const RealVector&), 
                                     void *fio, bool (*fi)(void*, const RealVector&), 
                                     std::vector<Curve> &ext_curve);

        virtual void extension_curve(const std::vector<RealVector> &curve, 
                                     const std::vector<RealVector> &domain_polygon, const std::vector<RealVector> &image_polygon, 
                                     std::vector<RealVector> &ext_curve);

        virtual void extension_curve(const Curve &curve, 
                                     const std::vector<RealVector> &domain_polygon,
                                     const std::vector<RealVector> &image_polygon,
                                     Curve &ext_curve);
};

#endif // _EXTENSION_


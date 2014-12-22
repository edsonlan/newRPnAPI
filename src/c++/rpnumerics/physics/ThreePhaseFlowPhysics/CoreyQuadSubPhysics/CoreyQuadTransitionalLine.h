#ifndef _COREYQUADTRANSITIONALLINE_
#define _COREYQUADTRANSITIONALLINE_

#include "RealVector.h"
#include "BifurcationCurve.h"

class CoreyQuadSubPhysics;

#define COREYQUADTRANSITIONALLINE_NONE (-1)
#define COREYQUADTRANSITIONALLINE_W_E    0
#define COREYQUADTRANSITIONALLINE_O_B    1
#define COREYQUADTRANSITIONALLINE_G_D    2

class CoreyQuadTransitionalLine: public BifurcationCurve {
    private:
    protected:
        CoreyQuadSubPhysics *corey_;

        static double WE_line(void *obj, const RealVector &p);
        static double OB_line(void *obj, const RealVector &p);
        static double GD_line(void *obj, const RealVector &p);
    public:
        CoreyQuadTransitionalLine(CoreyQuadSubPhysics *c);
        virtual ~CoreyQuadTransitionalLine();

        virtual void list_of_secondary_bifurcation_curves(std::vector<int> &type, 
                                                          std::vector<std::string> &name,
                                                          std::vector<void*> &object,
                                                          std::vector<double (*)(void*, const RealVector &)> &function){
            type.clear();
            name.clear();
            object.clear();
            function.clear();

            type.push_back(COREYQUADTRANSITIONALLINE_NONE);
            name.push_back(std::string("None"));
            object.push_back(0);
            function.push_back(0);

            type.push_back(COREYQUADTRANSITIONALLINE_W_E);
            name.push_back(std::string("W-E"));
            object.push_back(this);
            function.push_back(&WE_line);

            type.push_back(COREYQUADTRANSITIONALLINE_O_B);
            name.push_back(std::string("O-B"));
            object.push_back(this);
            function.push_back(&OB_line);

            type.push_back(COREYQUADTRANSITIONALLINE_G_D);
            name.push_back(std::string("G-D"));
            object.push_back(this);
            function.push_back(&GD_line);

            return;
        }

        void curve(int type, std::vector<RealVector> &c);
};

#endif // _COREYQUADTRANSITIONALLINE_


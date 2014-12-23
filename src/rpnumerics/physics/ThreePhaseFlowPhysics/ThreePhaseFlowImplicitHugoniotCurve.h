#ifndef _THREEPHASEFLOWIMPLICITHUGONIOTCURVE_
#define _THREEPHASEFLOWIMPLICITHUGONIOTCURVE_

#include "ImplicitHugoniotCurve.h"
#include "ThreePhaseFlowSubPhysics.h"

#define THREEPHASEFLOWIMPLICITHUGONIOTCURVE_G_VERTEX      1
#define THREEPHASEFLOWIMPLICITHUGONIOTCURVE_W_VERTEX      2
#define THREEPHASEFLOWIMPLICITHUGONIOTCURVE_O_VERTEX      3
#define THREEPHASEFLOWIMPLICITHUGONIOTCURVE_GW_SIDE       4
#define THREEPHASEFLOWIMPLICITHUGONIOTCURVE_WO_SIDE       5
#define THREEPHASEFLOWIMPLICITHUGONIOTCURVE_GO_SIDE       6

class ThreePhaseFlowImplicitHugoniotCurve : public ImplicitHugoniotCurve {
    private:
    protected:
        ThreePhaseFlowSubPhysics *subphysics_;
        ThreePhaseFlowPermeability *permeability_;

        double vel, muw, muo, /*mug,*/ grw, gro, grg;

        double rho_w_o, rho_o_w;
        double rho_w_g, rho_g_w;
        double rho_o_g, rho_g_o;

        double lambda_g_ref, lambda_o_ref, lambda_w_ref;

        static double gas_vertex_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state);
        static double water_vertex_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state);
        static double oil_vertex_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state);

        static double water_gas_side_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state);
        static double water_oil_side_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state);
        static double oil_gas_side_function(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state);

        double (*implicit_Hugoniot_function)(ThreePhaseFlowImplicitHugoniotCurve *obj, const RealVector &state);

        static int classical_function_on_square(ImplicitHugoniotCurve *obj, double *foncub, int i, int j);
        static int threephase_function_on_square(ImplicitHugoniotCurve *obj, double *foncub, int i, int j);

        int (*fonsq)(ImplicitHugoniotCurve *obj, double *foncub, int i, int j);
    public:
        ThreePhaseFlowImplicitHugoniotCurve(ThreePhaseFlowSubPhysics *s);
        virtual ~ThreePhaseFlowImplicitHugoniotCurve();

        int function_on_square(double *foncub, int i, int j);

        void curve(const ReferencePoint &ref, int type, std::vector<Curve> &c);

        void list_of_reference_points(std::vector<int> &type, std::vector<std::string> &name) const {
            type.clear();
            name.clear();

            type.push_back(IMPLICITHUGONIOTCURVE_GENERIC_POINT);
            name.push_back(std::string("Generic point (base class)"));

            type.push_back(THREEPHASEFLOWIMPLICITHUGONIOTCURVE_G_VERTEX);
            name.push_back(std::string("G vertex"));

            type.push_back(THREEPHASEFLOWIMPLICITHUGONIOTCURVE_W_VERTEX);
            name.push_back(std::string("W vertex"));

            type.push_back(THREEPHASEFLOWIMPLICITHUGONIOTCURVE_O_VERTEX);
            name.push_back(std::string("O vertex"));

            type.push_back(THREEPHASEFLOWIMPLICITHUGONIOTCURVE_GW_SIDE);
            name.push_back(std::string("GW side"));

            type.push_back(THREEPHASEFLOWIMPLICITHUGONIOTCURVE_WO_SIDE);
            name.push_back(std::string("WO side"));

            type.push_back(THREEPHASEFLOWIMPLICITHUGONIOTCURVE_GO_SIDE);
            name.push_back(std::string("GO side"));

            return;
        }

};

#endif // _THREEPHASEFLOWIMPLICITHUGONIOTCURVE_


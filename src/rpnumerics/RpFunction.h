/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RpFunction.h
 **/

#ifndef RPFUNCTION_H
#define RPFUNCTION_H


#define RP_SUCCESSFUL 1

/*
 ** ---------------------------------------------------------------
 ** Includes:
 */
//#include "ReturnCodes.h"
#include "JetMatrix.h"
#include "WaveState.h"
#include "Parameter.h"

/*
 ** ---------------------------------------------------------------
 ** Definitions:
 */


/*!@brief  The RpFunction class defines a generic function prototype.
 *
 *
 * f:u C Rm -> v C Rn
 *
 * The jet method returns the nth derivative
 * 
 * TODO:
 * make the jet method return f by default. In order to do that we must have
 * RealVector v parameter for f as a static or dynamic cast !!!
 * 
 *
 * @ingroup rpnumerics
 */
class RpFunction {
    protected:
        int Ceqs, vars; // Number of conserved variables and equations.

        std::vector<Parameter*> parameter;
    public:
        virtual ~RpFunction();
//        virtual RpFunction * clone() const = 0;
    
        /*! m coordinates function evaluation at u
         *this is the nth derivative calculation that might be available or not
         */

        virtual int jet(const WaveState &u, JetMatrix &m, int degree) const = 0;
        virtual int jet(const WaveState &u, JetMatrix &m_full, JetMatrix &m_redux, int degree) const {
            int info = jet(u, m_full, degree);

            m_redux = m_full;

            return info;
        }

        int jet(const RealVector &u, JetMatrix &m, int degree){
            return jet(WaveState(u), m, degree);
        }

        virtual void fill_with_jet(int n, const double *in, int degree, double *F, double *J, double *H) const;

        virtual int number_of_conservation_equations() const {return Ceqs;};
        virtual int number_of_variables() const {return vars;};

        virtual void parameters(std::vector<Parameter*> &vp){
            vp = parameter;
            return;
        }
};



#endif //RPFUNCTION_H

#ifndef _RIEMANNPROBLEM_
#define _RIEMANNPROBLEM_

#include "WaveCurve.h"

class RiemannProblem {
    private:
    protected:
        void half_profile(const WaveCurve &wavecurve, int curve, int point, int family,
                          std::vector<RealVector> &phase_state, std::vector<double> &speed);
    public:
        RiemannProblem();
        ~RiemannProblem();

        void profile(const WaveCurve &wavecurve1, int curve1, int point1, int family1,
                     const WaveCurve &wavecurve2, int curve2, int point2, int family2,
                     std::vector<RealVector> &phase_state, std::vector<double> &speed);

        void all_increase_profile(const WaveCurve &wavecurve1, int curve1, int point1, int family1,
                                          const WaveCurve &wavecurve2, int curve2, int point2, int family2,
                                          std::vector<RealVector> &phase_state, std::vector<double> &speed);
};

#endif // _RIEMANNPROBLEM_


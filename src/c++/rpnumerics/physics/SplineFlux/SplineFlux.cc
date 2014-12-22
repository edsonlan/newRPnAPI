#include "SplineFlux.h"

SplineFlux::SplineFlux(const std::vector<RealVector> &p) : FluxFunction() {
    int n = p.size();
    int n_base_func = n;

    for (int i = 0; i < n; i++) points.push_back(p[i]);

    ap::real_1d_array x, y;
    x.setlength(n);
    y.setlength(n);

    double xtemp, ytemp;

    for (int i = 0; i < n; i++){
        x(i) = points[i](0);
        y(i) = points[i](1);
    }

    // Generate the spline
    int spline_info;
    spline1dfitreport report;
    spline1dfitcubic(x, y, n, n_base_func,
                     spline_info, 
                     spline, 
                     report);
}

SplineFlux::~SplineFlux(){
}

SplineFlux * SplineFlux::clone() const {
    return new SplineFlux(points);
}

int SplineFlux::jet(const WaveState &u, JetMatrix &m, int degree) const {
    double large = 1000.0;


    double x = u(0);
    double y = large*u(1);

    double h, d_h, d2_h;
    spline1ddiff(spline, x, h, d_h, d2_h);

    m.resize(2, 2);

    if (degree >= 0){
        m.set(0, h);
        m.set(1, y);

        if (degree >= 1){
            m.set(0, 0, d_h);
            m.set(0, 1, 0.0);
            m.set(1, 0, 0.0);
            m.set(1, 1, large);

            if (degree >= 2){
                m.set(0, 0, 0, d2_h);
                m.set(0, 0, 1, 0.0);
                m.set(0, 1, 0, 0.0);
                m.set(0, 1, 1, 0.0);

                m.set(1, 0, 0, 0.0);
                m.set(1, 0, 1, 0.0);
                m.set(1, 1, 0, 0.0);
                m.set(1, 1, 1, 0.0);
            }
        }
    }

    return 2;
}

